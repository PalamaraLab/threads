# This file is part of the Threads software suite.
# Copyright (C) 2024 Threads Developers.
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
import h5py
import logging
import pgenlib
import importlib

os.environ["RAY_DEDUP_LOGS"] = "0"
import ray
import numpy as np

from threads_arg import (
    ThreadsLowMem,
    Matcher,
    ViterbiPath
)
from .utils import (
    make_recombination_from_map_and_pgen,
    make_constant_recombination_from_pgen,
    split_list,
    parse_demography,
    parse_region_string
)
from datetime import datetime

logger = logging.getLogger(__name__)


def serialize_paths(paths, positions, out, start=None, end=None):
    """
    Compress a list of paths across input positions
    """
    samples = [i for i in range(len(paths))]
    num_threads = len(samples)
    num_sites = len(positions)

    region_start = positions[0] if start == None else max(positions[0], start)
    region_end = positions[-1] + 1 if end == None else min(positions[-1] + 1, end + 1)

    thread_starts = []
    mut_starts = []
    thread_start = 0
    mut_start = 0
    all_bps, all_ids, all_ages, all_het_sites = [], [], [], []

    for i, path in enumerate(paths):
        path.map_positions(positions.astype(int))
        bps, ids, ages = path.dump_data_in_range(region_start, region_end)
        het_site_indices = np.array(path.het_sites, dtype=int)
        het_sites = positions[het_site_indices]
        het_indices_out = list(het_site_indices[(region_start <= het_sites) * (het_sites < region_end)])
        if i > 0:
            try:
                assert bps[0] == region_start
                if start is not None:
                    assert bps[0] >= start
                assert len(np.unique(bps)) == len(bps)
                assert bps[-1] <= positions[-1]
            except AssertionError:
                logger.info("Error serializing a segment with the following data:")
                logger.info(f"bps[0] {bps[0]}")
                logger.info(f"region_start {region_start}")
                logger.info(f"len(np.unique(bps)) {len(np.unique(bps))}")
                logger.info(f"len(bps) {len(bps)}")
                logger.info(f"bps[-1] {bps[-1]}")
                logger.info(f"positions[-1] {positions[-1]}")

        all_bps += bps
        all_ids += ids
        all_ages += ages
        all_het_sites += het_indices_out
        thread_starts.append(thread_start)
        mut_starts.append(mut_start)
        thread_start += len(bps)
        mut_start += len(het_indices_out)

    num_stitches = len(all_bps)
    num_mutations = len(all_het_sites)

    assert len(all_bps) == len(all_ids) == len(all_ages)
    assert len(paths) == len(samples)

    f = h5py.File(out, "w")
    f.attrs['datetime_created'] = datetime.now().isoformat()

    compression_opts = 9
    dset_samples = f.create_dataset("samples", (num_threads, 3), dtype=int, compression='gzip',
                                    compression_opts=compression_opts)
    dset_pos = f.create_dataset("positions", (num_sites), dtype=int, compression='gzip',
                                    compression_opts=compression_opts)
    # First L columns are random samples for imputation
    dset_targets = f.create_dataset("thread_targets", (num_stitches, 2), dtype=int, compression='gzip',
                                    compression_opts=compression_opts)
    dset_ages = f.create_dataset("thread_ages", (num_stitches), dtype=np.double, compression='gzip',
                                    compression_opts=compression_opts)
    dset_het_s = f.create_dataset("het_sites", (num_mutations), dtype=int, compression='gzip',
                                    compression_opts=compression_opts)
    dset_range = f.create_dataset("arg_range", (2), dtype=np.double, compression='gzip',
                                  compression_opts=compression_opts)

    dset_samples[:, 0] = samples
    dset_samples[:, 1] = thread_starts
    dset_samples[:, 2] = mut_starts
    dset_pos[:] = positions

    dset_targets[:, 0] = all_ids

    dset_targets[:, 1] = all_bps

    dset_ages[:] = all_ages

    dset_het_s[:] = all_het_sites

    dset_range[:] = [region_start, region_end]

    f.close()


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
        for b in range(n_batches):
            b_start = b * BATCH_SIZE
            b_end = min(M, (b+1) * BATCH_SIZE)
            g_size = b_end - b_start
            alleles_out = np.empty((g_size, num_samples_hap), dtype=np.int32)
            phased_out = np.empty((g_size, num_samples_hap // 2), dtype=np.uint8)
            reader.read_alleles_and_phasepresent_range(b_start, b_end, alleles_out, phased_out)
            for g in alleles_out:
                TLM.process_site_hets(g)

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
def threads_infer(pgen, map, recombination_rate, demography, mutation_rate, query_interval, match_group_interval, mode, num_threads, region, max_sample_batch_size, out):
    """Infer an ARG from genotype data"""
    start_time = time.time()
    logger.info(f"Starting Threads-infer with the following parameters:")
    logger.info(f"  pgen:                  {pgen}")
    logger.info(f"  map:                   {map}")
    logger.info(f"  recombination_rate:    {recombination_rate}")
    logger.info(f"  region:                {region}")
    logger.info(f"  demography:            {demography}")
    logger.info(f"  mutation_rate:         {mutation_rate}")
    logger.info(f"  query_interval:        {query_interval}")
    logger.info(f"  match_group_interval:  {match_group_interval}")
    logger.info(f"  num_threads:           {num_threads}")
    logger.info(f"  max_sample_batch_size: {max_sample_batch_size}")
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

    reader = pgenlib.PgenReader(pgen.encode())
    num_samples = reader.get_raw_sample_ct()
    num_sites = reader.get_variant_ct()
    logger.info(f"Will build an ARG on {2 * num_samples} haplotypes using {num_sites} sites")
    if max_sample_batch_size is None:
        max_sample_batch_size = 2 * num_samples

    # Initialize read batching
    M = len(physical_positions)
    assert num_sites == M
    BATCH_SIZE = int(4e7 // num_samples)
    n_batches = int(np.ceil(M / BATCH_SIZE))

    logger.info("Finding singletons")
    # Get singleton filter for the matching step
    alleles_out = None
    phased_out = None
    ac_mask = np.zeros(genetic_positions.shape, dtype=bool)
    for b in range(n_batches):
        b_start = b * BATCH_SIZE
        b_end = min(M, (b+1) * BATCH_SIZE)
        g_size = b_end - b_start
        alleles_out = np.empty((g_size, 2 * num_samples), dtype=np.int32)
        phased_out = np.empty((g_size, num_samples), dtype=np.uint8)
        reader.read_alleles_and_phasepresent_range(b_start, b_end, alleles_out, phased_out)
        # Filter out unphased variants and singletons
        ac_mask[b_start:b_end] = (alleles_out.sum(axis=1) > 1) * ~np.any(phased_out == 0, axis=1)

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

    # Matching step
    for b in range(n_batches):
        b_start = b * BATCH_SIZE
        b_end = min(M, (b+1) * BATCH_SIZE)
        g_size = b_end - b_start
        batch_mask = ac_mask[b_start:b_end]
        alleles_out = np.empty((g_size, 2 * num_samples), dtype=np.int32)
        phased_out = np.empty((g_size, num_samples), dtype=np.uint8)
        reader.read_alleles_and_phasepresent_range(b_start, b_end, alleles_out, phased_out)
        for g in alleles_out[batch_mask]:
            matcher.process_site(g)
    # Add top matches from adjacent sites to each match-chunk
    matcher.propagate_adjacent_matches()

    # From here we parallelise if we can
    actual_num_threads = min(len(os.sched_getaffinity(0)), num_threads)
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
    logger.info(f"Writing to {out}")
    # Save results
    serialize_paths(paths, physical_positions.astype(int), out, start=out_start, end=out_end)
    logger.info(f"Done in (s): {time.time() - start_time}")
