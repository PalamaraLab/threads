import logging
import multiprocessing
import os
import numpy as np
import pandas as pd

from tqdm import tqdm
from cyvcf2 import VCF
from threads_arg import ThreadsFastLS, ImputationMatcher
from scipy.sparse import csr_array, lil_matrix
from datetime import datetime
from typing import Dict, Tuple, List
from dataclasses import dataclass

from .fwbw import fwbw
from .utils import timer_block, log_nth_element

logger = logging.getLogger(__name__)


# Get available processes from os (rather than all with multiprocessing.cpu_count)
PROCESS_COUNT = len(os.sched_getaffinity(0))

@dataclass
class RecordMemo:
    id: int
    genotypes: None
    pos: None
    ref: None
    alt: None
    af: float

RecordMemoDict = Dict[str, RecordMemo]


class WriterVCF:
    """
    Custom VCF writer for imputation

    FIXME Review with Arni, use this customised version or alternative library
    """
    def __init__(self, filename):
        self.file = open(filename, "w")

    def write_header(self, samples, contig):
        f = self.file
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

    def write_site(self, genotypes, record, imputed, contig):
        imp_str = "IMP;" if imputed else ""
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

        f = self.file
        f.write(("\t".join([contig, pos, snp_id, ref, alt, qual, filter, f"{imp_str}AF={af:.4f}", "GT:DS", "\t".join(gt_strings)]) + "\n"))


def read_map_gz(map_gz):
    """
    Reading in map file for Li-Stephens

    Note that this is in the shapeit5 format not like the other maps we've been using
    """
    maps = pd.read_table(map_gz, sep='\s+')
    cm_pos = maps.cM.values.astype(np.float64)
    phys_pos = maps.pos.values.astype(np.float64)
    for i in range(1, len(cm_pos)):
        if cm_pos[i] <= cm_pos[i-1]:
            cm_pos[i] = cm_pos[i-1] + 1e-5
    return phys_pos, cm_pos


def parse_demography(demography):
    d = pd.read_table(demography, sep='\s+', header=None)
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


def site_arg_probability(site_posterior, imputation_thread, mutation_mapping, carriers, pos):
    """
    For each target sequence, the LS-fwbw finds out which of the panel sequences
    are most closely related to the target, rather than weighing by posterior.

    This is described under "Threading-based imputation" in the Methods section
    of the paper.
    """
    arg_probability = np.ones(site_posterior.shape)
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
            mut_lower, mut_upper = mutation_mapping.get_boundaries(s_id)
            mut_height = (mut_upper + mut_lower) / 2
            lam = 2. / height
            arg_probability[s_id] = 1 - np.exp(-lam * mut_height) * (1 + lam * mut_height)

    return arg_probability


class MutationMap:
    def __init__(self, snp, flipped, mapping_str):
        self.boundaries = {}
        self.snp = snp
        self.flipped = flipped == 1
        self.mapped = mapping_str != "NaN"
        self.uniquely_mapped = self.mapped and len(mapping_str.split(";")) == 1
        if self.mapped:
            for mut in mapping_str.split(";"): # FIXME expensive approx 10%
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

    def is_mapped(self, snp):
        try:
            return self.mut_dict[snp].is_mapped()
        except KeyError:
            return False

    def get_mapping(self, snp):
        return self.mut_dict[snp]


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


def _memoize_nth_record_process(filename, region, proc_idx, proc_max) -> List[Tuple[int, RecordMemo]]:
    """
    Process for for reading a region for a filename for every nth record

    Enumerating the enture VCF region per-process may seem counterintuitive
    because it's iterating many more times. However, the iteration itself is not
    expensive, it's the access to record.genotypes and conversion to a numpy
    array. Each process only calls these methods by filtering every record index
    against the process index.

    The returned list of memos is wrapped in a tuple that contains the original
    record index. This allows the final code to reconstruct the memos in order.

    proc_idx and proc_max are which process index this is and the max number of
    processes respectively.
    """
    results = []
    vcf = VCF(filename)
    for i, record in enumerate(vcf(region)):
        # This process only work on every nth record offset by i
        if i % proc_max == proc_idx:
            # The genotypes accessor and the conversion of the result (a python
            # list) into a flat bool numpy are the most expensive parts of this
            # operation, hence splitting into separate processes.
            genotypes = record.genotypes
            genotypes_flat = np.array(genotypes, dtype=bool)[:, :2].flatten()

            # Store all data required for imputation in a memo
            af = int(record.INFO["AC"]) / int(record.INFO["AN"])
            memo = RecordMemo(
                record.ID,
                genotypes_flat,
                record.POS,
                record.REF,
                record.ALT,
                af
            )

            # To avoid race conditions, the results are not injected directly
            # into a map, instead added to list with index to be sorted later
            # and store in map in original order.
            results.append((i, memo))

    return results


def _memoize_vcf_region_records(filename, region, cpu_count=PROCESS_COUNT) -> RecordMemoDict:
    """
    Given a VCF filename and region, generate a dictionary of record memos
    containing just the data required for imputation.
    """
    # Split the expensive parts reading records over available processes.
    # 'imemos' is shorthand for indexed memos, a list of tuples with index and memo
    shortname = os.path.basename(filename)
    with timer_block(f"memoising VCF {shortname}, region {region} ({cpu_count} CPUs)"):
        jobs_args = [(filename, region, i, cpu_count) for i in range(cpu_count)]
        with multiprocessing.Pool(processes=cpu_count) as pool:
            imemos = pool.starmap(_memoize_nth_record_process, jobs_args)

    # Flatten results (list of lists of tuples) into list of tuples
    imemos_flattened = [tup for tups in imemos for tup in tups]

    # Sort results by record index, the first element in tuple
    imemos_sorted = sorted(imemos_flattened, key=lambda tup: tup[0])

    # Collate into a dict looking up by memo id, using results_sorted ensure
    # that the dict are stored in the order the records were read from the file.
    record_dict = {memo.id: memo for _, memo in imemos_sorted}

    logger.info(f"Stored {len(record_dict)} memos from {shortname}")
    return record_dict


class Impute:
    """
    FIXME
    """
    def __init__(
        self,
        panel: str,
        target: str,
        map: str,
        mut: str,
        demography: str,
        out: str,
        region: str,
        mutation_rate=1.4e-8
    ):
        self.panel_dict = _memoize_vcf_region_records(panel, region)
        self.target_dict = _memoize_vcf_region_records(target, region)

        with timer_block("imputation"):
            target_samples = VCF(target).samples
            target_contigs = VCF(target).seqnames
            assert len(target_contigs) == 1
            chrom_num = target_contigs[0]

            # 2N as this is a collection of N diploid individuals, i.e. each has two haplotypes
            n_target_haps = 2 * len(target_samples)

            with timer_block("compute posteriors"):
                posteriors, imputation_threads = self._sparse_posteriors(
                    target=target,
                    map=map,
                    demography=demography,
                    region=region,
                    mutation_rate=mutation_rate
                )

            # FIXME work in progress
            # Cache special cases of first and last posteriors rows and an empty
            # table for all others that gets propaged when needed.
            row_len = len(posteriors[0].toarray())
            cached_posteriors_row_arrays = []
            cached_posteriors_first = []
            cached_posteriors_last = []
            for i, p in enumerate(posteriors):
                cached_posteriors_first.append(p[[0],:].toarray().flatten())
                cached_posteriors_last.append(p[[-1],:].toarray().flatten())
                cached_posteriors_row_arrays.append(None)

            def ensure_posteriors_cached(target_idx: int, snp_idx: int):
                if cached_posteriors_row_arrays[target_idx] is None:
                    cached_posteriors_row_arrays[target_idx] = [None] * row_len

                if cached_posteriors_row_arrays[target_idx][snp_idx] is None:
                    next_snp_row = posteriors[target_idx][[snp_idx],:].toarray()
                    cached_posteriors_row_arrays[target_idx][snp_idx] = (next_snp_row / np.sum(next_snp_row)).flatten()

            logger.info("Computing snps")
            phys_pos_array, _ = interpolate_map(target, map, region)
            num_snps = len(phys_pos_array)

            logger.info("Parsing mutations")
            mutation_container = MutationContainer(mut)

            # this will be the memory-heavy bit
            logger.info(f"Writing VCF header in {out}")
            vcf_writer = WriterVCF(out)
            vcf_writer.write_header(target_samples, chrom_num)

            snp_positions = []
            snp_ids = []
            with open(out, "a") as outfile:
                for record in VCF(target)(region):
                    snp_positions.append(record.POS)
                    snp_ids.append(record.ID)
                next_snp_idx = 0

                panel_vcf = VCF(panel)
                with timer_block(f"processing records"):
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
                            # FIXME conversion to float required for np.round
                            genotypes = np.array(self.target_snps[var_idx], dtype=float)
                        else:
                            while next_snp_idx < len(snp_positions) and pos >= snp_positions[next_snp_idx]:
                                next_snp_idx += 1
                                # FIXME Work in progress
                                # Rolling window of cached posteriors; they are computed when needed and dropped
                                # when out of range
                                for target_idx in range(n_target_haps):
                                    if not cached_posteriors_row_arrays[target_idx] is None:
                                        cache_window_end_idx = next_snp_idx - 2
                                        if cache_window_end_idx >= 0:
                                            if not cached_posteriors_row_arrays[target_idx][cache_window_end_idx] is None:
                                                cached_posteriors_row_arrays[target_idx][cache_window_end_idx] = None
                                if next_snp_idx == len(snp_positions):
                                    break

                            mutation_mapping = None
                            if mutation_container.is_mapped(var_id):
                                mutation_mapping = mutation_container.get_mapping(var_id)

                            genotypes = []
                            # FIXME WIP move out of loop?
                            # FIXME record.genotypes is not an ndarray, so conversion required. Custom version?
                            #panel_genotypes = np.array(record.genotypes, dtype=bool)[:, :2].flatten()
                            panel_genotypes = self.panel_dict[var_id].genotypes # FIXME remove dict lookup when switched to record list
                            carriers = (1 - panel_genotypes).nonzero()[0] if flipped else panel_genotypes.nonzero()[0]
                            for target_idx in range(n_target_haps):
                                # Extract and interpolate posterior
                                site_posterior = None
                                if next_snp_idx == 0:
                                    site_posterior = cached_posteriors_first[target_idx]
                                elif next_snp_idx == num_snps:
                                    site_posterior = cached_posteriors_last[target_idx]
                                else:
                                    ensure_posteriors_cached(target_idx, next_snp_idx)
                                    ensure_posteriors_cached(target_idx, next_snp_idx - 1)

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

                                    cached_target_idx = cached_posteriors_row_arrays[target_idx]
                                    # FIXME WIP move out of loop, and is np.average better in bulk?
                                    # - e.g. site_posterior = np.average([prev_target, next_target], axis=0, weights=[prev_wt, next_wt])
                                    prev_target = cached_target_idx[next_snp_idx - 1]
                                    next_target = cached_target_idx[next_snp_idx]
                                    site_posterior = prev_wt * prev_target + next_wt * next_target

                                # FIXME WIP block processing?
                                if mutation_mapping:
                                    arg_mask = site_arg_probability(site_posterior, imputation_threads[target_idx], mutation_mapping, carriers, pos)
                                    site_posterior = site_posterior * arg_mask

                                genotypes.append(np.sum(site_posterior, where=panel_genotypes))

                        genotypes = np.round(genotypes, decimals=3)
                        assert 0 <= np.max(genotypes) <= 1

                        vcf_writer.write_site(genotypes, record, imputed, chrom_num)


    def _sparse_posteriors(self, target, map, demography, region, mutation_rate):
        with timer_block("collating snps"):
            target_variants = set(self.target_dict.keys())
            self.panel_snps = np.array([record.genotypes for record in self.panel_dict.values() if record.id in target_variants])
            self.target_snps = np.array([record.genotypes for record in self.target_dict.values()])
            assert len(self.panel_snps) == len(self.target_snps)

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

        logger.info("Building panel...")
        for i, h in enumerate(self.panel_snps.transpose()):
            log_nth_element("Inserting haplotype", i)
            bwt.insert(h)

        with timer_block("reference matching"):
            ref_matches = reference_matching(self.panel_snps, self.target_snps, cm_pos_array)
            num_samples_panel = self.panel_snps.shape[1]

        mutation_rate = 0.0001
        cm_sizes = list(cm_pos_array[1:] - cm_pos_array[:-1])
        cm_sizes = np.array(cm_sizes + [cm_sizes[-1]])
        Ne = 20_000
        recombination_rates = 1 - np.exp(-4 * Ne * 0.01 * cm_sizes / num_samples_panel)

        posteriors = []
        imputation_threads = []
        L = 16
        logger.info("Doing posteriors...")
        for i, h_target in enumerate(self.target_snps.transpose()): # FIXME keep tdqm?
            with timer_block(f"impute and viterbi {i + 1}"):
                # Imputation thread with divergence matching
                imputation_thread = bwt.impute(list(h_target), L)
                imputation_threads.append(imputation_thread)
                matched_samples_viterbi = set([match_id for seg in imputation_thread for match_id in seg.ids])

                # All locally sampled matches
                matched_samples_matcher = (ref_matches[num_samples_panel + i])

                # Union of the two match sets
                matched_samples = np.array(list(matched_samples_viterbi.union(matched_samples_matcher)))

            with timer_block(f"fwbw {i + 1}"):
                posterior = fwbw(self.panel_snps[:, matched_samples], h_target[None, :], recombination_rates, mutation_rate)
                sparse_posterior = sparsify_posterior(posterior, matched_samples, num_snps=num_snps, num_samples=num_samples_panel)

                assert sparse_posterior.shape == (num_snps, num_samples_panel)
                posteriors.append(sparse_posterior)

        return posteriors, imputation_threads
