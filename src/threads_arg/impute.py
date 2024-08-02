import click
import logging
import time
import os
# import pgenlib
from tqdm import tqdm
from cyvcf2 import VCF
import h5py
from threads_arg import ThreadsFastLS, ImputationMatcher
import numpy as np
import pandas as pd
import lshmm
from scipy.sparse import csr_array, lil_matrix
from datetime import datetime

# this is only for imputation
def write_vcf_header(out, samples, contig):
    with open(out, "w") as f:
        f.write("##fileformat=VCFv4.2\n")
        f.write("##FILTER=<ID=PASS,Description=\"All filters passed\">\n")
        f.write(f"##fileDate={datetime.now().strftime('%m/%d/%Y, %H:%M:%S')}\n")
        f.write("##source=Threads 0.0\n")
        f.write(f"##contig=<ID={contig}>\n")
        f.write("##FPLOIDY=2\n")
        f.write("##INFO=<ID=IMP,Number=0,Type=Flag,Description=\"Imputed marker\">\n")
        f.write("##INFO=<ID=AF,Number=A,Type=Float,Description=\"Estimated allele frequency\">\n")
        f.write("##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Phased genotypes\">\n")
        f.write("##FORMAT=<ID=DS,Number=A,Type=Float,Description=\"Genotype dosage\">\n")
        f.write(("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t" + "\t".join(samples) + "\n"))

def write_site(f, genotypes, record, imputed, contig):
        imp_str = "IMP;" if imputed else ""
        # for (pos, snp_id, af, hap1_var, hap2_var, dosage_var) in zip(positions, snp_ids, afs, hap1_genotypes, hap2_genotypes, dosages):
        # imp_str = "" if pos in genotyped_pos else "IMP;"
        haps1 = genotypes[::2]
        haps2 = genotypes[1::2]
        dosages = haps1 + haps2
        pos = str(record.POS)
        snp_id = record.ID
        ref = record.REF
        alt = record.ALT
        assert len(alt) == 1
        alt = alt[0]
        qual = "."
        filter = "PASS"
        af = int(record.INFO["AC"]) / int(record.INFO["AN"])
        gt_strings = [f"{np.round(hap_1):.0f}|{np.round(hap_2):.0f}:{dosage:.3f}".rstrip("0").rstrip(".") for hap_1, hap_2, dosage in zip(haps1, haps2, dosages)]
        f.write(("\t".join([contig, pos, snp_id, ref, alt, qual, filter, f"{imp_str}AF={af:.4f}", "GT:DS", "\t".join(gt_strings)]) + "\n"))
        # f.write(("\t".join([contig, str(pos), snp_id.replace("'", "").replace(":", "_"), "A", "C", ".", "PASS", f"{imp_str}AF={af:.4f}", "GT:DS", "\t".join(gt_strings)]) + "\n"))

# def write_vcf(out_vcf_gz, samples, positions, snp_ids, afs, hap1_genotypes, hap2_genotypes, dosages, genotyped_pos):
#     chrom = "1"
#     # with gzip.open(out_vcf_gz, "wb") as f:
#     with open(out_vcf_gz, "w") as f:
#         f.write("##fileformat=VCFv4.2\n")
#         f.write("##FILTER=<ID=PASS,Description=\"All filters passed\">\n")
#         f.write(f"##fileDate={datetime.now().strftime('%m/%d/%Y, %H:%M:%S')}\n")
#         f.write("##source=Threads 0.0\n")
#         f.write(f"##contig=<ID={chrom}>\n")
#         f.write("##FPLOIDY=2\n")
#         f.write("##INFO=<ID=IMP,Number=0,Type=Flag,Description=\"Imputed marker\">\n")
#         f.write("##INFO=<ID=AF,Number=A,Type=Float,Description=\"Estimated allele frequency\">\n")
#         f.write("##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Phased genotypes\">\n")
#         f.write("##FORMAT=<ID=DS,Number=A,Type=Float,Description=\"Genotype dosage\">\n")
#         f.write(("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t" + "\t".join(samples) + "\n"))
#         for (pos, snp_id, af, hap1_var, hap2_var, dosage_var) in zip(positions, snp_ids, afs, hap1_genotypes, hap2_genotypes, dosages):
#             imp_str = "" if pos in genotyped_pos else "IMP;"
#             gt_strings = [f"{np.round(hap_1):.0f}|{np.round(hap_2):.0f}:{dosage:.3f}".rstrip("0").rstrip(".") for hap_1, hap_2, dosage in zip(hap1_var, hap2_var, dosage_var)]
#             f.write(("\t".join([chrom, str(pos), snp_id.replace("'", "").replace(":", "_"), "A", "C", ".", "PASS", f"{imp_str}AF={af:.4f}", "GT:DS", "\t".join(gt_strings)]) + "\n"))

logging.basicConfig(
    format='%(asctime)s %(levelname)-8s %(message)s',
    level=logging.INFO,
    datefmt='%Y-%m-%d %H:%M:%S')

def read_map_gz(map_gz):
    """
    Reading in map file for Li-Stephens
    THIS IS IN THE SHAPEIT5 FORMAT NOT LIKE THE OTHER MAPS WE'VE BEEN USING
    """
    # if (map_gz[:-3] == ".gz") :
    #     maps = pd.read_table(map_gz, header=None, compression='gzip')
    # else:
    #     maps = pd.read_table(map_gz, header=None)
    maps = pd.read_table(map_gz, delim_whitespace=True)
    
    cm_pos = maps.cM.values.astype(np.float64)
    phys_pos = maps.pos.values.astype(np.float64)
    for i in range(1, len(cm_pos)):
        if cm_pos[i] <= cm_pos[i-1]:
            cm_pos[i] = cm_pos[i-1] + 1e-5
    return phys_pos, cm_pos

def parse_demography(demography):
    d = pd.read_table(demography, delim_whitespace=True, header=None)
    return list(d[0]), list(d[1])

def interpolate_posterior(posterior, array_coords, wgs_coords):
    n_samples = posterior.shape[1]
    interpolated = np.empty((n_samples, len(wgs_coords)), dtype=float)
    interpolated[:] = np.nan
    for i in range(n_samples):
        interpolated[i, :] = np.interp(wgs_coords, array_coords, posterior[:, i])
    assert np.isnan(interpolated).sum() == 0

    return interpolated

def reference_matching(haps_panel, haps_target, cm_pos):
    num_reference = haps_panel.shape[1]
    num_target = haps_target.shape[1]
    matcher = ImputationMatcher(num_reference, num_target, cm_pos, 0.02, 4)
    all_genotypes = np.concatenate([haps_panel, haps_target], axis=1)
    for g in all_genotypes:
        matcher.process_site(g)
    return matcher.get_matches()

# def read_panel_wgs(wgs_path, matched_samples, start_idx, end_idx):
#     # funky mapping haploid diploid and back
#     diploid_ids = np.unique(matched_samples // 2)
#     all_haploid_ids = np.sort(np.concatenate((2 * diploid_ids, 2 * diploid_ids + 1)))
#     map_back = [np.searchsorted(all_haploid_ids, s_hap) for s_hap in matched_samples]

#     alleles_out = np.empty((end_idx - start_idx, len(all_haploid_ids)), dtype=np.int32)
#     phased_out = np.empty((end_idx - start_idx, len(diploid_ids)), dtype=np.uint8)

#     reader = pgenlib.PgenReader(wgs_path.encode(), sample_subset=diploid_ids.astype(np.uint32))
#     reader.read_alleles_and_phasepresent_range(start_idx, end_idx, alleles_out, phased_out)

#     return alleles_out.transpose()[map_back]

def site_arg_mask(site_posterior, imputation_thread, mutation_mapping, carriers, pos):
    arg_mask = np.ones(site_posterior.shape)
    num_segs = len(imputation_thread)
    segment = None
    if pos < imputation_thread[0].seg_start:
        segment = imputation_thread[0]
    else:
        seg_idx = 0
        while seg_idx < num_segs - 1 and imputation_thread[seg_idx + 1].seg_start < pos:
            seg_idx += 1
        segment = imputation_thread[seg_idx]

    for s_id, height in zip(segment.ids, segment.ages):
        if site_posterior[s_id] > 0 and s_id in carriers:
            try:
                mut_lower, mut_upper = mutation_mapping.get_boundaries(s_id)
            except KeyError:
                import pdb
                pdb.set_trace()
            mut_height = (mut_upper + mut_lower) / 2
            lam = 2. / height
            arg_mask[s_id] = 1 - np.exp(-lam * mut_height) * (1 + lam * mut_height)
    return arg_mask

        # if panel_wgs[panel_idx_map[s_id], i] == 1:
        #     mut_lower, mut_upper = mutation.get_boundaries(s_id)
        #     mut_height = (mut_upper + mut_lower) / 2
        #     lam = 2. / height

    # seg.ages
    # seg.ids
    # seg.seg_start
    # seg.weights


def build_arg_mask(imputation_thread, panel_wgs, phys_pos_out, mutation_container, panel_idx_map, start, end, arg_threshold=0.02):
    arg_mask = np.ones(panel_wgs.shape, dtype=float)
    allele_freqs = panel_wgs.mean(axis=0)
    current_seg_idx = 0
    for i, (pos, freq) in enumerate(zip(phys_pos_out, allele_freqs)):
        if mutation_container.mutation_mapped_at(pos):
            mutation = mutation_container.mutation_at(pos)

            if (not mutation.flipped and freq > arg_threshold) or (mutation.flipped and freq < 1 - arg_threshold):
                continue
            while current_seg_idx < len(imputation_thread) - 1 and pos >= imputation_thread[current_seg_idx + 1].seg_start:
                current_seg_idx += 1
            segment = imputation_thread[current_seg_idx]
            seg_start = max(start, segment.seg_start)
            seg_end = end if current_seg_idx == len(imputation_thread) - 1 else min(end, imputation_thread[current_seg_idx + 1].seg_start)
            try:
                assert seg_start <= pos < seg_end or pos == end
            except AssertionError:
                import pdb
                pdb.set_trace()

            # segment_height = segment.ages[0]
            try:
                focal_genotypes = panel_wgs[[panel_idx_map[s] for s in segment.ids], i]
            except KeyError:
                import pdb
                pdb.set_trace()
            if not mutation.flipped and focal_genotypes.sum() > 0:
                # something to check
                for s_id, height in zip(segment.ids, segment.ages):
                    if panel_wgs[panel_idx_map[s_id], i] == 1:
                        mut_lower, mut_upper = mutation.get_boundaries(s_id)
                        mut_height = (mut_upper + mut_lower) / 2
                        lam = 2. / height

                        # cdf_upper = np.exp(-lam * mut_upper) * (1 + lam * mut_upper)
                        # cdf_lower = np.exp(-lam * mut_lower) * (1 + lam * mut_lower)
                        # mut_distr = np.linspace(mut_lower, mut_upper, 10)

                        # dosage_lower = 0
                        # dosage_edge = (cdf_upper - cdf_lower) * (1 - np.mean(np.exp(-lam * mut_distr) * (1 + lam * mut_distr)))
                        # dosage_upper = 1 - cdf_upper

                        arg_mask[panel_idx_map[s_id], i] = 1 - np.exp(-lam * mut_height) * (1 + lam * mut_height)

                        # arg_mask[panel_idx_map[s_id], i] = 1 - np.exp(-lam * height) * (1 + lam * height)
                        # compute dosage here:
                        # if height >= upper:
                        #     continue
                        # elif height <= lower:
                        #     arg_mask[panel_idx_map[s_id], i] = 0
                        # else:
                            # arg_mask[panel_idx_map[s_id], i] = (height - lower) / (upper - lower)
                        # arg_mask[panel_idx_map[s], i]
            elif mutation.flipped and focal_genotypes.sum() < len(focal_genotypes):
                # something to check
                for s_id, height in zip(segment.ids, segment.ages):
                    if panel_wgs[[panel_idx_map[s_id]], i] == 0:
                        mut_lower, mut_upper = mutation.get_boundaries(s_id)
                        mut_height = (mut_upper + mut_lower) / 2
                        lam = 2. / height

                        arg_mask[panel_idx_map[s_id], i] = 1 - np.exp(-lam * mut_height) * (1 + lam * mut_height)

    return arg_mask

class MutationMap:
    def __init__(self, snp, flipped, mapping_str):
        self.boundaries = {}
        self.snp = snp
        self.flipped = flipped == 1
        self.mapped = mapping_str != "NaN"
        self.uniquely_mapped = self.mapped and len(mapping_str.split(";")) == 1
        if self.mapped:
            for mut in mapping_str.split(";"):
                ids, lower, upper = mut.split(",")
                for sample_id in ids.split("."):
                    self.boundaries[int(sample_id)] = (float(lower), float(upper))

    def is_carrier(self, sample_id):
        return int(sample_id) in self.boundaries

    def is_mapped(self):
        return self.mapped

    def get_boundaries(self, sample_id):
        """If uniquely mapped (with -1), assume we're querying for a carrier, otherwise look up mapping"""
        if self.uniquely_mapped:
            return self.boundaries[-1]
        else:
            return self.boundaries[int(sample_id)]

class MutationContainer:
    def __init__(self, mut_path):
        self.mut_dict = {}
        with open(mut_path, 'r') as mutfile:
            for line in mutfile:
                snp, pos, flipped, mapping_str = line.strip().split()
                self.mut_dict[snp] = MutationMap(snp, int(flipped), mapping_str)
                # self.mut_dict[int(pos)] = MutationMap(snp, int(flipped), mapping_str)

    def is_mapped(self, snp):
        try:
            return self.mut_dict[snp].is_mapped()
        except KeyError:
            return False

    def get_mapping(self, snp):
        return self.mut_dict[snp]
    
    # def mutation_mapped_at(self, snp):
    #     try:
    #         # return self.mut_dict[int(pos)].mapped
    #         return self.mut_dict[snp].mapped
    #     except KeyError:
    #         return False

def get_variants(target, region):
    vcf = VCF(target)
    variants = []
    for record in vcf(region):
        variants.append(record.ID)
    return variants

def read_panel_snps(panel, target, region):
    # read genotypes from panel that are in target
    target_variants = get_variants(target, region)
    # panel_variants = get_variants(panel, region)

    genotypes = []
    vcf = VCF(panel)
    for record in vcf(region):
        if record.ID in target_variants:
            genotypes.append(np.array(record.genotypes)[:, :2].flatten())
    return np.array(genotypes)

def read_target_snps(target, region):
    genotypes = []
    vcf = VCF(target)
    for record in vcf(region):
        genotypes.append(np.array(record.genotypes)[:, :2].flatten())
    return np.array(genotypes)

def sparsify_posterior(P, matched_samples, num_snps, num_samples):
    """
    Expand to a compressed n_snps x n_samples matrix
    """
    assert P.shape == (num_snps, len(matched_samples))
    matrix = lil_matrix(np.zeros(shape=(num_snps, num_samples)))
    P[P <= 1 / num_samples] = 0
    for i, p in enumerate(P):
        assert np.sum(p) > 0
        q = p / np.sum(p)
        for j in np.nonzero(q)[0]:
            matrix[i, matched_samples[j]] = q[j]
    return csr_array(matrix)

def interpolate_map(vcf_path, map_path, region):
    # get vcf positions
    vcf = VCF(vcf_path)
    positions = []
    for record in vcf(region):
        positions.append(record.POS)
    positions = np.array(positions)

    map_bp, map_cm = read_map_gz(map_path)
    cm_array = np.interp(positions, map_bp, map_cm)
    return positions, cm_array

def sparse_posteriors(panel, target, map, demography, region, mutation_rate):
    panel_snps = read_panel_snps(panel, target, region)
    target_snps = read_target_snps(target, region)
    assert len(panel_snps) == len(target_snps)

    phys_pos_array, cm_pos_array = interpolate_map(target, map, region)
    num_snps = len(phys_pos_array)

    sparse_sites = True
    use_hmm = False
    ne_times, ne_sizes = parse_demography(demography)
    bwt = ThreadsFastLS(phys_pos_array,
                        cm_pos_array,
                        mutation_rate,
                        ne_sizes,
                        ne_times,
                        sparse_sites,
                        use_hmm=use_hmm)

    logging.info("Building panel...")
    for h in panel_snps.transpose():
        bwt.insert(h)

    ref_matches = reference_matching(panel_snps, target_snps, cm_pos_array)
    num_samples_panel = panel_snps.shape[1]
    num_samples_target = target_snps.shape[1]

    # Todo: get expected time from input demography :-)
    mutation_rates = np.array([0.0001] * num_snps)
    cm_sizes = list(cm_pos_array[1:] - cm_pos_array[:-1])
    cm_sizes = np.array(cm_sizes + [cm_sizes[-1]])
    Ne = 20_000
    recombination_rates = 1 - np.exp(-4 * Ne * 0.01 * cm_sizes / num_samples_panel)

    posteriors = []  # he, he
    imputation_threads = []  # he, he
    L = 16
    logging.info("Doing posteriors")
    for i, h_target in tqdm(enumerate(target_snps.transpose())):
        # Imputation thread with divergence matching
        imputation_thread = bwt.impute(list(h_target), L)
        imputation_threads.append(imputation_thread)
        # imputation_thread = bwt.impute(h_target, L=16)
        matched_samples_viterbi = set([match_id for seg in imputation_thread for match_id in seg.ids])

        # All locally sampled matches
        matched_samples_matcher = (ref_matches[num_samples_panel + i])

        # Union of the two match sets
        matched_samples = np.array(list(matched_samples_viterbi.union(matched_samples_matcher)))
        panel_idx_map = {s: i for i, s in enumerate(matched_samples)}

        # Janky hmm (should re-do this more efficiently, this won't scale)
        forward_array, normalisation_factor_from_forward, log_likelihood = lshmm.forwards(panel_snps[:, matched_samples], h_target[None, :], recombination_rates, None, mutation_rates, scale_mutation_based_on_n_alleles=False)
        backward_array = lshmm.backwards(panel_snps[:, matched_samples], h_target[None, :], normalisation_factor_from_forward, recombination_rates, None, mutation_rates, scale_mutation_based_on_n_alleles=False)
        posterior = forward_array * backward_array

        sparse_posterior = sparsify_posterior(posterior, matched_samples, num_snps=num_snps, num_samples=num_samples_panel)

        assert sparse_posterior.shape == (num_snps, num_samples_panel)
        posteriors.append(sparse_posterior)
    return posteriors, imputation_threads

# @click.command()
# @click.option("--panel_array", required=True, help="pgen array panel. Array!!")
# @click.option("--panel_wgs", required=True, help="pgen array panel. Wgs!!")
# @click.option("--target", required=True, help="pgen array targets")
# @click.option("--mutation_file", required=True, help="pgen array targets", default=None)
# @click.option("--map_gz_array", required=True, help="Path to genotype map")
# @click.option("--map_gz_wgs", required=True, help="Path to genotype map")
# @click.option("--mutation_rate", required=True, type=float, help="Per-site-per-generation SNP mutation rate.")
# @click.option("--demography", required=True, help="Path to file containing demographic history.")
# @click.option("--out", required=True, help="Path to output .vcf file.")
# @click.option("--out_arg", required=True, help="Path to output .vcf file.")
# @click.option("--l", default=128, help="Maximum number of matching states sampled per segment. Default is 128 (too high? too high)")
# @click.option("--start", type=int)
# @click.option("--end", type=int)
# def impute(panel_array, panel_wgs, target, mutation_file, map_gz_array, map_gz_wgs, mutation_rate, demography, out, out_arg, l, start, end):
@click.command()
@click.option("--panel", required=True, help="pgen array panel. Array!!")
@click.option("--target", required=True, help="pgen array targets")
@click.option("--mut", required=True, help="pgen array targets")
@click.option("--map", required=True, help="Path to genotype map")
@click.option("--mutation_rate", type=float, help="Per-site-per-generation SNP mutation rate.", default=1.4e-8)
@click.option("--demography", required=True, help="Path to file containing demographic history.")
@click.option("--out", required=True, help="Path to output .vcf file.")
@click.option("--region", required=True, type=str)
def impute(panel, target, map, mut, demography, out, region, mutation_rate=1.4e-8):
    threads_impute(panel, target, map, mut, demography, out, region, mutation_rate)


def threads_impute(panel, target, map, mut, demography, out, region, mutation_rate=1.4e-8):
    start_time = time.time()

    panel_samples = VCF(panel).samples
    target_samples = VCF(target).samples
    target_contigs = VCF(target).seqnames
    assert len(target_contigs) == 1
    chrom_num = target_contigs[0]
    n_panel_haps = 2 * len(panel_samples)
    n_target_haps = 2 * len(target_samples)

    logging.info("Starting!")
    posteriors, imputation_threads = sparse_posteriors(panel=panel,
        target=target,
        map=map,
        demography=demography,
        region=region,
        mutation_rate=mutation_rate)
    logging.info("Done with posteriors!")

    phys_pos_array, cm_pos_array = interpolate_map(target, map, region)
    num_snps = len(phys_pos_array)

    logging.info("Parsing mutations!")
    mutation_container = MutationContainer(mut)
    logging.info("Done with mutations!")

    # this will be the memory-heavy bit
    logging.info(f"Imputing and writing to {out}")
    write_vcf_header(out, target_samples, chrom_num)

    snp_positions = []
    snp_ids = []
    with open(out, "a") as outfile:
        for record in VCF(target)(region):
            snp_positions.append(record.POS)
            snp_ids.append(record.ID)
        target_genotypes = read_target_snps(target, region)
        next_snp_idx = 0

        panel_vcf = VCF(panel)
        for record in panel_vcf(region):
            var_id = record.ID
            pos = record.POS
            af = int(record.INFO["AC"]) / int(record.INFO["AN"])
            flipped = af > 0.5
            genotypes = None
            imputed = True
            if var_id in snp_ids:
                imputed = False
                var_idx = snp_ids.index(var_id)
                genotypes = target_genotypes[var_idx]
                # write genotype from there
            else:
                while next_snp_idx < len(snp_positions) and pos >= snp_positions[next_snp_idx]:
                    next_snp_idx += 1
                    if next_snp_idx == len(snp_positions):
                        break
                genotypes = []
                panel_genotypes = np.array(record.genotypes)[:, :2].flatten()
                carriers = (1 - panel_genotypes).nonzero()[0] if flipped else panel_genotypes.nonzero()[0]
                for target_idx in range(n_target_haps):
                    # Extract and interpolate posterior
                    site_posterior = None
                    if next_snp_idx == 0:
                        # FIXME replace _getrow with [] operator after profiling
                        site_posterior = posteriors[target_idx]._getrow(0).toarray()
                    elif next_snp_idx == num_snps:
                        site_posterior = posteriors[target_idx]._getrow(-1).toarray()
                    else:
                        bp_prev = snp_positions[next_snp_idx - 1]
                        bp_next = snp_positions[next_snp_idx]
                        assert bp_prev <= pos <= bp_next
                        prev_wt, next_wt = None, None
                        if bp_prev == bp_next:
                            prev_wt = 1
                            next_wt = 0
                        else:
                            prev_wt = (pos - bp_prev) / (bp_next - bp_prev)
                            next_wt = 1 - prev_wt
                        site_posterior = prev_wt * posteriors[target_idx]._getrow(next_snp_idx - 1).toarray() + next_wt * posteriors[target_idx]._getrow(next_snp_idx).toarray()
                        site_posterior /= np.sum(site_posterior)
                    site_posterior = site_posterior.flatten()

                    if mutation_container.is_mapped(var_id):
                        mutation_mapping = mutation_container.get_mapping(var_id)
                        arg_mask = site_arg_mask(site_posterior, imputation_threads[target_idx], mutation_mapping, carriers, pos)
                        genotypes.append((site_posterior * arg_mask * panel_genotypes).sum())
                    else:
                        genotypes.append((site_posterior * panel_genotypes).sum())
            genotypes = np.round(genotypes, decimals=3)
            assert 0 <= np.max(genotypes) <= 1
            write_site(outfile, genotypes, record, imputed, chrom_num)
    logging.info(f"Done! in {time.time() - start_time:.3f} (s)")
    # posterior = interpolate_posterior(posterior, cm_pos_array, cm_pos_out)
    # panel_alleles = read_panel_wgs(panel_wgs, matched_samples, start_idx, end_idx)
    # assert posterior.shape == panel_alleles.shape
    # arg_mask = build_arg_mask(imputation_thread, panel_alleles, phys_pos_out, mutation_container, panel_idx_map, start=start, end=end)
    # assert arg_mask.shape == panel_alleles.shape

    # imputed_dosages[i] = (posterior * panel_alleles).sum(axis=0)
    # arg_imputed_dosages[i] = (posterior * arg_mask * panel_alleles).sum(axis=0)
    
    # assert np.isnan(imputed_dosages).sum() == 0
    # assert np.isnan(arg_imputed_dosages).sum() == 0

    # snp_ids = [f"SNP_{p}" for p in phys_pos_out]
    # afs = [0 for p in phys_pos_out]
    # samples = [f"Sample_{p}" for p in range(num_samples_target // 2)]

    # hap1_genotypes = imputed_dosages[::2]
    # hap2_genotypes = imputed_dosages[1::2]
    # dosages =  (imputed_dosages[::2] + imputed_dosages[1::2]).transpose()
    # hap1_genotypes = np.rint(imputed_dosages[::2]).astype(int).transpose()
    # hap2_genotypes = np.rint(imputed_dosages[1::2]).astype(int).transpose()
    # assert len(phys_pos_out) == len(hap1_genotypes) == len(hap2_genotypes) == len(dosages) == len(snp_ids)
    # write_vcf(out, samples, phys_pos_out, snp_ids, afs, hap1_genotypes, hap2_genotypes, dosages, phys_pos_array)

    # arg_hap1_genotypes = arg_imputed_dosages[::2]
    # arg_hap2_genotypes = arg_imputed_dosages[1::2]
    # arg_dosages =  (arg_imputed_dosages[::2] + arg_imputed_dosages[1::2]).transpose()
    # arg_hap1_genotypes = np.rint(arg_imputed_dosages[::2]).astype(int).transpose()
    # arg_hap2_genotypes = np.rint(arg_imputed_dosages[1::2]).astype(int).transpose()
    # assert len(phys_pos_out) == len(arg_hap1_genotypes) == len(arg_hap2_genotypes) == len(arg_dosages) == len(snp_ids)
    # write_vcf(out_arg, samples, phys_pos_out, snp_ids, afs, arg_hap1_genotypes, arg_hap2_genotypes, arg_dosages, phys_pos_array)

if __name__ == "__main__":
    impute()
