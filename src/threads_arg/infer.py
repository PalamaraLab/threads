# This file is part of the Threads software suite.
# Copyright (C) 2024-2025 Threads Developers.
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <http://www.gnu.org/licenses/>.

import os
import gc
import time
import logging
import pgenlib
import importlib

os.environ["RAY_DEDUP_LOGS"] = "0"
import ray
import numpy as np
import pandas as pd

from threads_arg import (
    ThreadsLowMem,
    Matcher,
    ViterbiPath,
    ThreadingInstructions,
    ConsistencyWrapper,
    GenotypeIterator
)
from .utils import (
    make_recombination_from_map_and_pgen,
    make_constant_recombination_from_pgen,
    split_list,
    parse_demography,
    iterate_pgen,
    read_positions_and_ids,
    parse_region_string,
    default_process_count,
    read_variant_metadata,
    read_sample_names
)

from .serialization import serialize_instructions

from .allele_ages import estimate_ages

logger = logging.getLogger(__name__)


def partial_viterbi(pgen, mode, num_samples_hap, physical_positions, genetic_positions, demography, mu, sample_batch, s_match_group, match_cm_positions, max_sample_batch_size, num_threads, thread_id):
    """Parallelized ARG inference sub-routine"""
    start_time = time.time()
    # Force logging to go straight to stdout, instead of into ray tmp files
    logging.shutdown()
    importlib.reload(logging)
    pid = os.getpid()
    logging.basicConfig(
        format=f"%(asctime)s %(levelname)-8s PID {pid} %(message)s",
        level=logging.INFO,
        datefmt='%Y-%m-%d %H:%M:%S'
    )
    local_logger = logging.getLogger(__name__)

    reader = pgenlib.PgenReader(pgen.encode())
    ne_times, ne = parse_demography(demography)

    sparse = None
    if mode == "array":
        sparse = True
    elif mode == "wgs":
        sparse = False
    else:
        raise RuntimeError

    # Batching here saves a small amount of memory 
    num_samples = len(sample_batch)
    sample_indices = list(range(num_samples))
    num_subsets = int(np.ceil(num_samples / max_sample_batch_size))
    sample_batch_subsets = split_list(sample_batch, num_subsets)
    sample_index_subsets = split_list(sample_indices, num_subsets)

    seg_starts = []
    match_ids  = []
    heights    = []
    hetsites   = []

    # Main loop, infer threading instructions for each sample in the batch
    for batch_index, (sample_batch_subset, sample_index_subset) in enumerate(zip(sample_batch_subsets, sample_index_subsets)):
        local_logger.info(f"Thread {thread_id}: HMM batch {batch_index + 1}/{num_subsets}...")
        TLM = ThreadsLowMem(sample_batch_subset, physical_positions, genetic_positions, ne, ne_times, mu, sparse)

        # Warning: this creates big copies of data
        if num_subsets == 1:
            TLM.initialize_viterbi(s_match_group, match_cm_positions)
        else:
            TLM.initialize_viterbi([[s[k] for k in sample_index_subset] for s in s_match_group], match_cm_positions)

        M = reader.get_variant_ct()
        BATCH_SIZE = int(4e7 // num_samples_hap)
        n_batches = int(np.ceil(M / BATCH_SIZE))

        # Initialize pruning parameters
        a_counter   = 0
        prune_threshold = 10 * num_samples_hap
        prune_count = 0
        last_prune  = 0

        # Iterate across the genotypes and run Li-Stephens inference
        for b in range(n_batches):
            # Read genotypes and check for phase
            b_start = b * BATCH_SIZE
            b_end = min(M, (b+1) * BATCH_SIZE)
            g_size = b_end - b_start
            alleles_out = np.empty((g_size, num_samples_hap), dtype=np.int32)
            phased_out = np.empty((g_size, num_samples_hap // 2), dtype=np.uint8)
            if (phased_out == 0).any():
                unphased_sites, unphased_samples = (1 - phased_out).nonzero()
                ALLELE_UNPHASED_HET = -7
                alleles_out[unphased_sites, 2 * unphased_samples] = ALLELE_UNPHASED_HET
                alleles_out[unphased_sites, 2 * unphased_samples + 1] = ALLELE_UNPHASED_HET
            reader.read_alleles_and_phasepresent_range(b_start, b_end, alleles_out, phased_out)

            # For each variant in chunk, pass the genotypes through Threads-LS
            for g in alleles_out:
                TLM.process_site_viterbi(g)

                # Regularly prune the number of open branches and reset the pruning threshold
                if a_counter % 10 == 0:
                    n_branches = TLM.count_branches()
                    if n_branches > prune_threshold:
                        if a_counter - last_prune <= 30:
                            prune_threshold *= 2
                        TLM.prune()
                        prune_count += 1
                        last_prune = a_counter
                a_counter += 1

        # Construct paths
        TLM.traceback()

        # Add heterozygous sites to each path segment
        iterate_pgen(pgen, lambda _, g: TLM.process_site_hets(g))

        # Add coalescence time to each segment
        TLM.date_segments()

        # Serialize paths and append to batch data
        seg_starts_sub, match_ids_sub, heights_sub, hetsites_sub = TLM.serialize_paths()
        seg_starts += seg_starts_sub
        match_ids += match_ids_sub
        heights += heights_sub
        hetsites += hetsites_sub

    end_time = time.time()
    local_logger.info(f"Thread {thread_id}: HMMs done in (s) {end_time - start_time:.1f}")
    return seg_starts, match_ids, heights, hetsites


# Implementation is separated from Click entrypoint for use in tests
def threads_infer(pgen, map, recombination_rate, demography, mutation_rate, fit_to_data, allele_ages, query_interval, match_group_interval, mode, num_threads, region, max_sample_batch_size, save_metadata, out):
    """Infer an ARG from genotype data"""
    start_time = time.time()
    logger.info(f"Starting Threads-infer with the following parameters:")
    logger.info(f"  pgen:                  {pgen}")
    logger.info(f"  map:                   {map}")
    logger.info(f"  recombination_rate:    {recombination_rate}")
    logger.info(f"  region:                {region}")
    logger.info(f"  demography:            {demography}")
    logger.info(f"  mutation_rate:         {mutation_rate}")
    logger.info(f"  allele_ages:           {allele_ages}")
    logger.info(f"  fit_to_data:           {fit_to_data}")
    logger.info(f"  query_interval:        {query_interval}")
    logger.info(f"  match_group_interval:  {match_group_interval}")
    logger.info(f"  num_threads:           {num_threads}")
    logger.info(f"  max_sample_batch_size: {max_sample_batch_size}")
    logger.info(f"  save_metadata:         {save_metadata}")
    logger.info(f"  out:                   {out}")

    out_start = None
    out_end = None
    chrom = None
    if region is not None:
        chrom, out_start, out_end = parse_region_string(region)

    # Initialize region, genetic map, and genotype reader
    if map is not None:
        logger.info(f"Using recombination rates from {map}")
        genetic_positions, physical_positions = make_recombination_from_map_and_pgen(map, pgen, chrom)
    else:
        logger.info(f"Using constant recombination rate of {recombination_rate}")
        genetic_positions, physical_positions = make_constant_recombination_from_pgen(pgen, recombination_rate, chrom)

    # Load/set CHR, POS, ID, REF, ALT, QUAL, FILTER
    variant_metadata = read_variant_metadata(pgen) if save_metadata else None

    # Load sample names
    sample_names = read_sample_names(pgen) if save_metadata else None

    if fit_to_data and (physical_positions[1:] - physical_positions[:-1] <= 0).any():
        raise RuntimeError("Sites must be strictly increasing when --fit-to-data is set.")


    reader = pgenlib.PgenReader(pgen.encode())
    num_samples = reader.get_raw_sample_ct()
    num_sites = reader.get_variant_ct()
    assert num_sites == len(physical_positions)
    logger.info(f"Will build an ARG on {2 * num_samples} haplotypes using {num_sites} sites")
    if max_sample_batch_size is None:
        max_sample_batch_size = 2 * num_samples

    logger.info("Finding singletons")
    # Get singleton filter for the matching step
    alleles_out = None
    phased_out = None
    ac_mask = []
    iterate_pgen(pgen, lambda i, g: ac_mask.append(1 < g.sum() < 2 * num_samples))
    ac_mask = np.array(ac_mask, dtype=bool)
    assert ac_mask.shape == genetic_positions.shape

    logger.info("Running PBWT matching")
    # There are four params here to mess with:
    #   - query interval
    #   - match group size
    #   - neighborhood size
    #   - min_matches
    # Keep min_matches low for small sample sizes, can be 2 up to ~1,000 but then > 3
    # Smaller numbers run faster and consume less memory
    MIN_MATCHES = 4
    neighborhood_size = 4
    matcher = Matcher(2 * num_samples, genetic_positions[ac_mask], query_interval, match_group_interval, neighborhood_size, MIN_MATCHES)
    def matcher_callback(i, g, mask, matcher):
        if mask[i]:
            matcher.process_site(g)
    iterate_pgen(pgen, matcher_callback, mask=ac_mask, matcher=matcher)

    # Add top matches from adjacent sites to each match-chunk
    matcher.propagate_adjacent_matches()

    # From here we parallelise if we can
    actual_num_threads = min(default_process_count(), num_threads)
    logger.info(f"Requested {num_threads} threads, found {actual_num_threads}.")
    paths = []
    if actual_num_threads > 1:
        # Warning: this creates big copies, these matches are the main source of memory usage
        sample_batches = split_list(list(range(2 * num_samples)), actual_num_threads)
        match_cm_positions = matcher.cm_positions()

        del alleles_out
        del phased_out
        gc.collect()
        partial_viterbi_remote = ray.remote(partial_viterbi)
        ray.init()
        # Parallelised threading instructions
        results = ray.get([partial_viterbi_remote.remote(
            pgen,
            mode,
            2 * num_samples,
            physical_positions,
            genetic_positions,
            demography,
            mutation_rate,
            sample_batch,
            matcher.serializable_matches(sample_batch),
            match_cm_positions,
            max_sample_batch_size,
            actual_num_threads,
            thread_id) for thread_id, sample_batch in enumerate(sample_batches)])
        ray.shutdown()
        # Combine results from each thread
        for sample_batch, result_tuple in zip(sample_batches, results):
            for sample_id, seg_starts, match_ids, heights, hetsites in zip(sample_batch, *result_tuple):
                paths.append(ViterbiPath(sample_id, seg_starts, match_ids, heights, hetsites))
    else:
        sample_batch = list(range(2 * num_samples))
        s_match_group = matcher.serializable_matches(sample_batch)
        match_cm_positions = matcher.cm_positions()
        matcher.clear()
        del matcher
        gc.collect()
        thread_id = 1
        # Single-threaded threading instructions
        results = partial_viterbi(
            pgen,
            mode,
            2 * num_samples,
            physical_positions,
            genetic_positions,
            demography,
            mutation_rate,
            sample_batch,
            s_match_group,
            match_cm_positions,
            max_sample_batch_size,
            actual_num_threads,
            thread_id)

        for sample_id, seg_starts, match_ids, heights, hetsites in zip(sample_batch, *results):
            paths.append(ViterbiPath(sample_id, seg_starts, match_ids, heights, hetsites))

    # Inference is done, now we slice up the threading instructions and keep only the region requested
    region_start = physical_positions[0] if out_start is None else max(physical_positions[0], out_start)
    region_end = physical_positions[-1] + 1 if out_end is None else min(physical_positions[-1] + 1, out_end + 1)
    instructions = ThreadingInstructions(paths, int(region_start), int(region_end), physical_positions.astype(int))

    if fit_to_data:
        logger.info("Starting data-consistency post-processing")

        allele_age_estimates = None
        if allele_ages is None:
            logger.info("Inferring allele ages from data")
            allele_age_estimates = estimate_ages(instructions, 3 * actual_num_threads, actual_num_threads)
            assert len(allele_age_estimates) == len(instructions.positions)
        else:
            logger.info(f"Reading allele ages from {allele_ages}")
            allele_age_estimates = []
            _, ids = read_positions_and_ids(pgen)
            age_table = pd.read_table(allele_ages, header=None, names=["SNP", "POS", "AGE"])
            age_table = age_table[age_table["SNP"].astype(str).isin(ids)]
            allele_age_estimates = age_table["AGE"].values
            try:
                assert age_table.shape[0] == len(instructions.positions) == len(allele_age_estimates)
            except AssertionError:
                raise RuntimeError(f"Allele age estimates do not match markers in the region requested, expected {len(instructions.positions)} age estimates.")

        # Start the consistifying
        logger.info("Post-processing threading instructions to fit to data")
        gt_it = GenotypeIterator(instructions)
        cw = ConsistencyWrapper(instructions, allele_age_estimates)
        while gt_it.has_next_genotype():
            g = np.array(gt_it.next_genotype())
            cw.process_site(g)

        consistent_instructions = cw.get_consistent_instructions()
        logger.info(f"Writing to {out}")
        serialize_instructions(consistent_instructions,
                               out,
                               variant_metadata=variant_metadata,
                               allele_ages=allele_age_estimates,
                               sample_names=sample_names)
    else:
        logger.info(f"Writing to {out}")
        serialize_instructions(instructions,
                               out,
                               variant_metadata=variant_metadata,
                               sample_names=sample_names)

    # Save results
    logger.info(f"Done in (s): {time.time() - start_time}")

