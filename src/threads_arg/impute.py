import logging
import multiprocessing
import os
import numpy as np
import pandas as pd
import sys

from tqdm import tqdm
from cyvcf2 import VCF
from threads_arg import ThreadsFastLS, ImputationMatcher
from scipy.sparse import csr_array, lil_matrix
from datetime import datetime
from typing import Dict, Tuple, List, Union
from dataclasses import dataclass

# Use bisect_left to optimise when available (key search only in python >= 3.10)
BISECT_LEFT_KEY_SEARCH = sys.version_info[:3] >= (3, 10)
if BISECT_LEFT_KEY_SEARCH:
    from bisect import bisect_left

from .fwbw import fwbw
from .utils import (
    timer_block,
    TimerTotal,
    read_map_file,
    parse_demography,
    default_process_count
)

logger = logging.getLogger(__name__)

@dataclass
class RecordMemo:
    """
    Memo of the record fields required for imputation. Cached at the start of
    the run to avoid repeatedly re-evaluating records and data like allele
    frequency and genotypes. The latter is particular expensive so is converted
    from a python list into a numpy bool array.
    """
    id: int
    genotypes: None
    pos: None
    ref: None
    alt: None
    af: float

RecordMemoDict = Dict[str, RecordMemo]


class WriterVCF:
    """
    Custom VCF writer for imputation. Writes to stdout if filename is None
    """
    def __init__(self, filename: Union[str, None]):
        if filename:
            self.file = open(filename, "w")
        else:
            self.file = sys.stdout

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
        pos = str(record.pos)
        snp_id = record.id
        ref = record.ref
        alt = record.alt
        af = record.af
        assert len(alt) == 1
        alt = alt[0]
        qual = "."
        filter = "PASS"
        gt_strings = [f"{np.round(hap_1):.0f}|{np.round(hap_2):.0f}:{dosage:.3f}".rstrip("0").rstrip(".") for hap_1, hap_2, dosage in zip(haps1, haps2, dosages)]
        f = self.file
        f.write(("\t".join([contig, pos, snp_id, ref, alt, qual, filter, f"{imp_str}AF={af:.4f}", "GT:DS", "\t".join(gt_strings)]) + "\n"))


def _reference_matching(haps_panel, haps_target, cm_pos):
    num_reference = haps_panel.shape[1]
    num_target = haps_target.shape[1]
    matcher = ImputationMatcher(num_reference, num_target, cm_pos, 0.02, 4)
    all_genotypes = np.concatenate([haps_panel, haps_target], axis=1)
    for g in tqdm(all_genotypes, mininterval=1):
        matcher.process_site(g)
    return matcher.get_matches()


def _active_site_arg_delta(
    active_site_posterior,
    active_indexes,
    imputation_thread,
    mutation_mapping,
    carriers,
    record
):
    """
    For each target sequence, the LS-fwbw finds out which of the panel sequences
    are most closely related to the target, rather than weighing by posterior.

    This is described under "Threading-based imputation" in the Methods section
    of the paper.
    """
    delta = 0

    # Find the nearest segment in imputation thread based on position
    if BISECT_LEFT_KEY_SEARCH:
        # Binary search for seg_idx before pos
        seg_idx = bisect_left(imputation_thread[1:], record.pos, key=lambda x: x.seg_start)
    else:
        # When key search not available, fall back to slower linear search
        num_segs = len(imputation_thread)
        seg_idx = 0
        while seg_idx < num_segs - 1 and imputation_thread[seg_idx + 1].seg_start < record.pos:
            seg_idx += 1
    segment = imputation_thread[seg_idx]

    for s_id, height in zip(segment.ids, segment.ages):
        # Reject anything that is not active (i.e. not True in record.genotypes)
        # and its index in active_site_posterior array
        if s_id in active_indexes:
            active_idx = active_indexes[s_id]

            if active_site_posterior[active_idx] > 0 and s_id in carriers:
                mut_lower, mut_upper = mutation_mapping.get_boundaries(s_id)
                mut_height = (mut_upper + mut_lower) / 2
                lam = 2. / height
                lam_mut = lam * mut_height
                arg_prob = 1 - np.exp(-lam_mut) * (1 + lam_mut)

                # Outside this function, the delta is applied to a pre-computed
                # sum of active site posteriors, so (arg_prob - 1) effectively
                # removes this sites contribution from that previously-summed
                # value in order to re-apply the probability-adjusted version.
                delta += active_site_posterior[active_idx] * (arg_prob - 1)

    return delta


class MutationMap:
    """
    Wrapper for the mapping of a variant to one or more branches in the arg,
    containing the descendants of the mutations and the boundaries of their
    mapped edges.
    """
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
    """
    A container for all MutationMap objects along the genomic chunk
    """
    def __init__(self, mut_path):
        self.mut_dict = {}
        with open(mut_path, 'r') as mutfile:
            for line in mutfile:
                snp, _pos, flipped, mapping_str = line.strip().split()
                self.mut_dict[snp] = MutationMap(snp, int(flipped), mapping_str)

    def is_mapped(self, snp):
        try:
            return self.mut_dict[snp].is_mapped()
        except KeyError:
            return False

    def get_mapping(self, snp):
        return self.mut_dict[snp]


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


def _memoize_nth_record_process_star(args):
    """
    Unpack star args for _memoize_nth_record_process

    This is required to get tqdm working with a multiprocess pool. The usual way
    would be with multiprocessing pool starmap(), but this does not work with
    tqdm. Instead imap() must be used with a shim method to unpack args.
    """
    return _memoize_nth_record_process(*args)


def _memoize_vcf_region_records(filename, region, process_count=None) -> RecordMemoDict:
    """
    Given a VCF filename and region, generate a dictionary of record memos
    containing just the data required for imputation.
    """
    if not process_count:
        process_count = default_process_count()

    # Split the expensive parts reading records over available processes.
    # 'imemos' is shorthand for indexed memos, a list of tuples with index and memo
    shortname = os.path.basename(filename)
    with timer_block(f"memoising VCF {shortname}, region {region} ({process_count} CPUs)", False):
        jobs_args = [(filename, region, i, process_count) for i in range(process_count)]
        with multiprocessing.Pool(processes=process_count) as pool:
            # To use tqdm with a pool, use imap with shim method to unpack args.
            imemos = list(tqdm(
                pool.imap(_memoize_nth_record_process_star, jobs_args),
                total=len(jobs_args),
                mininterval=1
            ))

    # Flatten results (list of sub lists of tuples) into list of tuples
    imemos_flattened = [
        tup
        for sub_list in imemos
        for tup in sub_list
    ]

    # Sort results by record index, the first element in tuple
    imemos_sorted = sorted(imemos_flattened, key=lambda tup: tup[0])

    # Collate into a dict looking up by memo id, using results_sorted ensure
    # that the dict are stored in the order the records were read from the file.
    record_dict = {memo.id: memo for _, memo in imemos_sorted}

    logger.info(f"Stored {len(record_dict)} memos from {shortname}")
    return record_dict


class CachedPosteriorSnps:
    """
    Maintain a row of posteriors per target index per snp index

    The main impute loop only deals with prev and next values so this only
    caches 2 by default.
    """
    def __init__(self, posteriors, max_size: int=2):
        self.posteriors = posteriors
        self.posteriors_by_snp_idx = {}
        self.recents = []
        self.max_size = max_size

    def __getitem__(self, snp_idx: int):
        """
        Index lookup will automatically cache and drop rows as required
        """
        # Flag snp index as being recently used
        self._flag_recent(snp_idx)

        # Reuse cached value if present
        if snp_idx in self.posteriors_by_snp_idx:
            return self.posteriors_by_snp_idx[snp_idx]

        # Rebuild snp data
        col_len = len(self.posteriors)
        row_len = self.posteriors[0].shape[1]
        target_posteriors = np.empty(shape=(col_len, row_len), dtype=np.float64)
        for i, p in enumerate(self.posteriors):
            posteriors = p[[snp_idx],:].toarray()
            target_posteriors[i] = posteriors / np.sum(posteriors)

        # Cache and clear out-of-date entries
        self.posteriors_by_snp_idx[snp_idx] = target_posteriors
        self._remove_oldest()

        return target_posteriors

    def _flag_recent(self, snp_idx):
        # Latest accessed snp index is moved to the front of the recents list
        if snp_idx in self.recents:
            self.recents.remove(snp_idx)
        self.recents.insert(0, snp_idx)

    def _remove_oldest(self):
        # Clear anything off the end of the recents list from the cache
        for idx in self.recents[self.max_size:]:
            del self.posteriors_by_snp_idx[idx]
        self.recents = self.recents[:self.max_size]


class Impute:
    """
    Wrapper class for the Threads imputation routine

    See paper for more details.
    """
    def __init__(
        self,
        panel: str,
        target: str,
        map: str,
        mut: str,
        demography: str,
        out: Union[str, None],
        region: str,
        mutation_rate=1.4e-8
    ):
        # If out is not specified then the generated file goes to stdout. These
        # changes ensure that progress bars and logger is not also printed.
        if not out:
            def disabled_tqdm(it, *_args, **_kvargs):
                return it
            global tqdm
            tqdm = disabled_tqdm
            logging.disable(logging.INFO)

        logger.info(f"Starting Threads-impute with parameters")
        logger.info(f"  panel:         {panel}")
        logger.info(f"  target:        {target}")
        logger.info(f"  map:           {map}")
        logger.info(f"  mut:           {mut}")
        logger.info(f"  demography:    {demography}")
        logger.info(f"  out:           {out or 'STDOUT'}")
        logger.info(f"  region:        {region}")
        logger.info(f"  mutation_rate: {mutation_rate}")

        with timer_block("impute"):
            self._load_records_and_snps(panel, target, map, region)

            self.target_samples = VCF(target).samples
            self.target_contigs = VCF(target).seqnames
            assert len(self.target_contigs) == 1
            self.chrom_num = self.target_contigs[0]

            # 2N as this is a collection of N diploid individuals, i.e. each has two haplotypes
            self.n_target_haps = 2 * len(self.target_samples)

            with timer_block("sparse posteriors"):
                self._sparse_posteriors(
                    demography=demography,
                    mutation_rate=mutation_rate
                )

            logger.info("Parsing mutations")
            self.mutation_container = MutationContainer(mut)

            logger.info("Init iteration")
            self._init_step_snp()
            self.posteriors_snp_cache = CachedPosteriorSnps(self.posteriors)

            logger.info(f"Writing VCF header to {out or 'STDOUT'}")
            self.vcf_writer = WriterVCF(out)
            self.vcf_writer.write_header(self.target_samples, self.chrom_num)

            self.tt_posteriors = TimerTotal("interpolate posteriors")
            self.tt_mutation_mapping = TimerTotal("mutation mapping")
            self.tt_genotypes = TimerTotal("genotypes")

            with timer_block(f"processing records"):
                self._process_and_write()

            logger.info(self.tt_posteriors)
            logger.info(self.tt_mutation_mapping)
            logger.info(self.tt_genotypes)


    def _load_records_and_snps(self, panel: str, target: str, map: str, region: str):
        """
        Called on initialisation to preload frequently-used values.
        """
        with timer_block("memozing VCF panel and region"):
            self.panel_dict = _memoize_vcf_region_records(panel, region)
            self.target_dict = _memoize_vcf_region_records(target, region)

        with timer_block("collating snps", False):
            target_variants = set(self.target_dict.keys())
            self.panel_snps = np.array([record.genotypes for record in self.panel_dict.values() if record.id in target_variants])
            self.target_snps = np.array([record.genotypes for record in self.target_dict.values()])
            assert len(self.panel_snps) == len(self.target_snps)
            self.num_samples_panel = self.panel_snps.shape[1]
            logger.info(f"Samples panel size {self.num_samples_panel}")

        with timer_block("getting VCF positions", False):
            self.phys_pos_array = np.array([record.pos for record in self.target_dict.values()])
            map_bp, map_cm, _chrom = read_map_file(map)
            self.cm_pos_array = np.interp(self.phys_pos_array, map_bp, map_cm)
            self.num_snps = len(self.phys_pos_array)
            logger.info(f"Pos array snps size {self.num_snps}")


    def _sparse_posteriors(self, demography, mutation_rate):
        """
        Compute sparsified Li-Stephens posterior matrices for all targets across
        the whole chunk
        """
        sparse_sites = True
        use_hmm = False
        ne_times, ne_sizes = parse_demography(demography)
        bwt = ThreadsFastLS(self.phys_pos_array,
                            self.cm_pos_array,
                            mutation_rate,
                            ne_sizes,
                            ne_times,
                            sparse_sites,
                            use_hmm=use_hmm)

        with timer_block("building panel", False):
            for h in tqdm(self.panel_snps.transpose(), mininterval=1):
                bwt.insert(h)

        with timer_block("reference matching"):
            ref_matches = _reference_matching(self.panel_snps, self.target_snps, self.cm_pos_array)

        mutation_rate = 0.0001
        cm_sizes = list(self.cm_pos_array[1:] - self.cm_pos_array[:-1])
        cm_sizes = np.array(cm_sizes + [cm_sizes[-1]])
        Ne = 20_000
        recombination_rates = 1 - np.exp(-4 * Ne * 0.01 * cm_sizes / self.num_samples_panel)

        posteriors = []
        imputation_threads = []
        L = 16
        with timer_block(f"posteriors"):
            tt_impute = TimerTotal("PBWT impute")
            tt_fwbw = TimerTotal("fwbw")
            target_transpose = self.target_snps.transpose()

            # tqdm total required to ensure bar updates correctly
            for i, h_target in tqdm(enumerate(target_transpose), total=len(target_transpose), mininterval=1):
                with tt_impute:
                    # Imputation thread with divergence matching
                    imputation_thread = bwt.impute(list(h_target), L)
                    imputation_threads.append(imputation_thread)
                    matched_samples_viterbi = set([match_id for seg in imputation_thread for match_id in seg.ids])

                    # All locally sampled matches
                    matched_samples_matcher = (ref_matches[self.num_samples_panel + i])

                    # Union of the two match sets
                    matched_samples = np.array(list(matched_samples_viterbi.union(matched_samples_matcher)))

                with tt_fwbw:
                    posterior = fwbw(self.panel_snps[:, matched_samples], h_target[None, :], recombination_rates, mutation_rate)
                    sparse_posterior = self._sparsify_posterior(posterior, matched_samples)
                    posteriors.append(sparse_posterior)

            logger.info(tt_impute)
            logger.info(tt_fwbw)

        self.posteriors = posteriors
        self.imputation_threads = imputation_threads


    def _sparsify_posterior(self, posterior, matched_samples):
        """
        Expand to a compressed n_snps x n_samples matrix
        """
        assert posterior.shape == (self.num_snps, len(matched_samples))
        matrix = lil_matrix((self.num_snps, self.num_samples_panel), dtype=np.float64)
        posterior[posterior <= 1 / self.num_samples_panel] = 0
        for i, p in enumerate(posterior):
            assert np.sum(p) > 0
            q = p / np.sum(p)
            for j in np.nonzero(q)[0]:
                matrix[i, matched_samples[j]] = q[j]
        return csr_array(matrix)


    def _init_step_snp(self):
        """
        Prepare any variables required for stepping through variants
        """
        self.snp_positions = [record.pos for record in self.target_dict.values()]
        self.snp_id_indexes = {record.id: i for i, record in enumerate(self.target_dict.values())}
        self.snp_idx = 0


    def _next_step_snp(self, record):
        """
        Interpolate and compute posterior weights for the current variant
        """
        site_posteriors = None
        mutation_mapping = None
        carriers = None

        pos = record.pos
        flipped = record.af > 0.5
        snp_idx = self.snp_idx

        def only_active(posteriors):
            return posteriors[:,record.genotypes]

        with self.tt_posteriors:
            # Scan for next snp index
            while snp_idx < self.num_snps and pos >= self.snp_positions[snp_idx]:
                snp_idx += 1
            self.snp_idx = snp_idx

            # Explicit values for first and last entries
            if snp_idx == 0:
                site_posteriors = only_active(self.posteriors_snp_cache[0])
            elif snp_idx == self.num_snps:
                site_posteriors = only_active(self.posteriors_snp_cache[-1])
            else:
                # Compute weights for position based on snp positions
                bp_prev = self.snp_positions[snp_idx - 1]
                bp_next = self.snp_positions[snp_idx]
                assert bp_prev <= pos <= bp_next
                prev_wt = (pos - bp_prev) / (bp_next - bp_prev) if bp_next != bp_prev else 1
                next_wt = 1 - prev_wt

                # Provide quick lookup to prev and next targets
                prev_target = only_active(self.posteriors_snp_cache[snp_idx - 1])
                next_target = only_active(self.posteriors_snp_cache[snp_idx])
                site_posteriors = prev_wt * prev_target + next_wt * next_target

            if self.mutation_container.is_mapped(record.id):
                mutation_mapping = self.mutation_container.get_mapping(record.id)
                carriers = (1 - record.genotypes).nonzero()[0] if flipped else record.genotypes.nonzero()[0]

            active_positions = np.where(record.genotypes)[0]
            active_indexes = {pos: i for i, pos in enumerate(active_positions)}

            return site_posteriors, mutation_mapping, carriers, active_indexes


    def _compute_genotypes(
        self,
        record,
        active_site_posteriors,
        active_indexes,
        mutation_mapping,
        carriers
    ):
        """
        Perform imputation for the current variant
        """
        def compute_delta(active_site_posterior, i):
            return _active_site_arg_delta(
                active_site_posterior,
                active_indexes,
                self.imputation_threads[i],
                mutation_mapping,
                carriers,
                record
            )

        genotypes = np.array([np.sum(asp) for asp in active_site_posteriors])
        if mutation_mapping:
            with self.tt_mutation_mapping:
                deltas = np.array([compute_delta(asp, i) for i, asp in enumerate(active_site_posteriors)])
                genotypes += deltas

        genotypes = np.round(genotypes, decimals=3)
        assert np.min(genotypes) >= 0
        assert np.max(genotypes) <= 1
        return genotypes


    def _process_and_write(self):
        """
        Process the loaded records and write resulting VCF
        """
        tt_write_vcf = TimerTotal("write vcf")

        for record in tqdm(self.panel_dict.values(), mininterval=1):
            genotypes = None
            imputed = True

            if record.id in self.snp_id_indexes:
                # If this variant is present on the genotyping array, then copy
                # it directly and note in the output genotypes that it was not
                # imputed
                imputed = False
                var_idx = self.snp_id_indexes[record.id]
                genotypes = np.array(self.target_snps[var_idx], dtype=float)
            else:
                # If this variant is not present on the genotyping array, then
                # run the Threads imputation routine
                active_site_posteriors, mutation_mapping, carriers, active_indexes = self._next_step_snp(record)

                genotypes = self._compute_genotypes(
                    record,
                    active_site_posteriors,
                    active_indexes,
                    mutation_mapping,
                    carriers,
                )

            with tt_write_vcf:
                self.vcf_writer.write_site(genotypes, record, imputed, self.chrom_num)

        logger.info(tt_write_vcf)
