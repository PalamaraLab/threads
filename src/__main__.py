import os
import gc
import sys
import time
import h5py
import tszip
import click
import logging
import pgenlib
import importlib
from memory_profiler import profile
import subprocess
import arg_needle_lib

os.environ["RAY_DEDUP_LOGS"] = "0"
import ray
import numpy as np
import xarray as xr
import pandas as pd

from threads import Threads, ThreadsLowMem, Matcher, ViterbiPath
from .utils import decompress_threads, interpolate_map, parse_demography, get_map_from_bim
from datetime import datetime
# from pandas_plink import read_plink1_bin
# import xarray as xr
from cyvcf2 import VCF, Writer

def print_help_msg(command):
    with click.Context(command) as ctx:
        click.echo(command.get_help(ctx))

logging.basicConfig(
    format='%(asctime)s %(levelname)-8s %(message)s',
    level=logging.INFO,
    datefmt='%Y-%m-%d %H:%M:%S')

def goodbye():
    # All credit to https://www.asciiart.eu/miscellaneous/dna
    print(
    """
    `-:-.   ,-;"`-:-.   ,-;"`-:-.   ,-;"`-:-.   ,-;"
       `=`,'=_______________________________=`,'=/
         y== ||Thank you for using Threads|| y==/  
       ,=,-<=‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾=,-<=`.
    ,-'-'   `-=_,-'-'   `-=_,-'-'   `-=_,-'-'   `-=_
    """
    )

def phase_distance(arg, arg_pos, target, het_carriers, hom_carriers):
    """This function only used for phasing. Compute 'phasing distance' as described in the paper."""
    if len(het_carriers) + 2 * len(hom_carriers) == 1:
        # negative so lower number (with longer pendant branch) wins
        return -arg.node(target).parent_edge_at(arg_pos).parent.height
    phase_distance = 0
    for carrier in het_carriers:
        phase_distance += min(arg.mrca(target, 2 * carrier, arg_pos).height, arg.mrca(target, 2 * carrier + 1, arg_pos).height)
    for hom in hom_carriers:
        phase_distance += arg.mrca(target, 2 * hom, arg_pos).height + arg.mrca(target, 2 * hom + 1, arg_pos).height
    return phase_distance

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
                logging.info("Error serializing a segment with the following data:")
                logging.info(f"bps[0] {bps[0]}")
                logging.info(f"region_start {region_start}")
                logging.info(f"len(np.unique(bps)) {len(np.unique(bps))}")
                logging.info(f"len(bps) {len(bps)}")
                logging.info(f"bps[-1] {bps[-1]}")
                logging.info(f"positions[-1] {positions[-1]}")

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

# def serialize(out, threads, samples, positions, arg_range, L=1):
#     num_threads = len(samples)
#     num_sites = len(positions)

#     thread_starts = []
#     mut_starts = []
#     thread_start = 0
#     mut_start = 0
#     all_bps, all_ids, all_ages, all_het_sites = [], [], [], []
#     for i, thread in enumerate(threads):
#         bps, ids, ages, het_sites = thread
#         all_bps += bps
#         all_ids += ids
#         all_ages += ages
#         all_het_sites += het_sites
#         thread_starts.append(thread_start)
#         mut_starts.append(mut_start)
#         thread_start += len(bps)
#         mut_start += len(het_sites)
#     num_stitches = len(all_bps)
#     num_mutations = len(all_het_sites)

#     assert len(all_bps) == len(all_ids) == len(all_ages)
#     assert len(threads) == len(samples)

#     f = h5py.File(out, "w")
#     f.attrs['datetime_created'] = datetime.now().isoformat()

#     compression_opts = 9
#     dset_samples = f.create_dataset("samples", (num_threads, 3), dtype=int, compression='gzip',
#                                     compression_opts=compression_opts)
#     dset_pos = f.create_dataset("positions", (num_sites), dtype=int, compression='gzip',
#                                  compression_opts=compression_opts)
#     # First L columns are random samples for imputation
#     dset_targets = f.create_dataset("thread_targets", (num_stitches, L + 2), dtype=int, compression='gzip',
#                                     compression_opts=compression_opts)
#     dset_ages = f.create_dataset("thread_ages", (num_stitches), dtype=np.double, compression='gzip',
#                                  compression_opts=compression_opts)
#     dset_het_s = f.create_dataset("het_sites", (num_mutations), dtype=int, compression='gzip',
#                                   compression_opts=compression_opts)
#     dset_range = f.create_dataset("arg_range", (2), dtype=np.double, compression='gzip',
#                                   compression_opts=compression_opts)

#     dset_samples[:, 0] = samples
#     dset_samples[:, 1] = thread_starts
#     dset_samples[:, 2] = mut_starts
#     dset_pos[:] = positions

#     for i in range(L + 1):
#         dset_targets[:, i] = [ids[i] for ids in all_ids]

#     dset_targets[:, L + 1] = all_bps

#     dset_ages[:] = all_ages

#     dset_het_s[:] = all_het_sites

#     dset_range[:] = arg_range

#     f.close()

def threads_to_arg(thread_dict, noise=0.0, max_n=None, verify=False):
    """
    Assemble threading instructions into an ARG
    """
    threading_instructions = thread_dict['threads']
    pos = thread_dict['positions']
    N = max_n if max_n is not None else np.max(thread_dict['samples'].astype(int)) + 1
    logging.info(f"Will thread {N} haploids")
    arg_start, arg_end = thread_dict['arg_range']
    # "+ 2" so we can include mutations on the last site, ARG is end-inclusive, we're not.
    arg = arg_needle_lib.ARG(0, arg_end - arg_start + 2, reserved_samples=N)
    arg.set_offset(int(arg_start))

    # How this should work:
    for i, t in enumerate(threading_instructions[:N]):
        if i % 1000 == 0:
            logging.info(f"Sequence {i + 1}...")
        arg.add_sample(str(i))
        if i > 0:
            section_starts, thread_ids, thread_heights = t
            if len(thread_ids.shape) == 2:
                thread_ids = thread_ids[:, 0]
            thread_heights += thread_heights * np.random.normal(0.0, noise, len(thread_heights))
            # this is weird, should fix
            arg_starts = [s - arg.offset for s in section_starts]
            try:
                if arg_starts[-1] >= arg.end:
                    arg.thread_sample([s - arg.offset for s in section_starts[:-1]], thread_ids[:-1], thread_heights[:-1])
                else:
                    arg.thread_sample([s - arg.offset for s in section_starts], thread_ids, thread_heights)
            except ValueError:
                import pdb
                pdb.set_trace()
    logging.info(f"Done threading")

    if verify:
        logging.info("Verifying ARG...")
        arg.check_basic()
    return arg

def split_list(list, n):
    """Yield n number of sequential chunks from l."""
    sublists = []
    d, r = divmod(len(list), n)
    for i in range(n):
        si = (d+1)*(i if i < r else r) + d*(0 if i < r else i - r)
        sublists.append(list[si:si+(d+1 if i < r else d)])
    return sublists

def partial_viterbi(pgen, mode, num_samples_hap, physical_positions, genetic_positions, demography, mu, sample_batch, s_match_group, match_cm_positions, max_sample_batch_size, num_threads, thread_id):
    """Parallelized ARG inference sub-routine"""
    start_time = time.time()
    # Force logging to go straight to stdout, instead of into ray tmp files
    logging.shutdown()
    importlib.reload(logging)
    pid = os.getpid()
    logging.basicConfig(format=f"%(asctime)s %(levelname)-8s PID {pid} %(message)s", 
                        level=logging.INFO,
                        datefmt='%Y-%m-%d %H:%M:%S')
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
                # We encode "unphased het" by "-7"
                alleles_out[unphased_sites, 2 * unphased_samples] = -7
                alleles_out[unphased_sites, 2 * unphased_samples + 1] = -7
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


@click.command()
@click.argument("command", type=click.Choice(["files", "infer", "convert", "impute"]))
@click.option("--pgen", required=True)
@click.option("--map_gz", default=None)
@click.option("--recombination_rate", default=1.3e-8, type=float)
@click.option("--demography", required=True)
@click.option("--mode", required=True, type=click.Choice(['array', 'wgs']))
@click.option("--query_interval", type=float, default=0.01, help="For the preliminary haplotype matching. Given in cM.")
@click.option("--match_group_interval", type=float, default=0.5, help="For the preliminary haplotype matching. Given in cM.")
@click.option("--mutation_rate", required=True, type=float)
@click.option("--num_threads", type=int, default=1)
@click.option("--region", help="Of format 123-345, end-exclusive") 
@click.option("--max_sample_batch_size", help="Of format 123-345, end-exclusive", default=None, type=int) 
@click.option("--out")
def infer(command, pgen, map_gz, recombination_rate, demography, mutation_rate, query_interval, match_group_interval, mode, num_threads, region, max_sample_batch_size, out):
    assert command == "infer"
    start_time = time.time()
    logging.info(f"Starting Threads-infer with the following parameters:")
    logging.info(f"  pgen:                  {pgen}")
    logging.info(f"  map_gz:                {map_gz}")
    logging.info(f"  recombination_rate:    {recombination_rate}")
    logging.info(f"  region:                {region}")
    logging.info(f"  demography:            {demography}")
    logging.info(f"  mutation_rate:         {mutation_rate}")
    logging.info(f"  query_interval:        {query_interval}")
    logging.info(f"  match_group_interval:  {match_group_interval}")
    logging.info(f"  num_threads:           {num_threads}")
    logging.info(f"  max_sample_batch_size: {max_sample_batch_size}")
    logging.info(f"  out:                   {out}")

    # Initialize region, genetic map, and genotype reader
    if map_gz is not None:
        logging.info(f"Using recombination rates from {map_gz}")
        genetic_positions, physical_positions = interpolate_map(map_gz, pgen)
    else:
        logging.info(f"Using constant recombination rate of {recombination_rate}")
        genetic_positions, physical_positions = get_map_from_bim(pgen, recombination_rate)

    out_start, out_end = None, None 
    if region is not None: 
        out_start, out_end = [int(r) for r in region.split("-")]

    reader = pgenlib.PgenReader(pgen.encode())
    num_samples = reader.get_raw_sample_ct()
    num_sites = reader.get_variant_ct()
    logging.info(f"Will build an ARG on {2 * num_samples} haplotypes using {num_sites} sites")
    if max_sample_batch_size is None:
        max_sample_batch_size = 2 * num_samples

    # Initialize read batching
    M = len(physical_positions)
    assert num_sites == M
    BATCH_SIZE = int(4e7 // num_samples)
    n_batches = int(np.ceil(M / BATCH_SIZE))

    logging.info("Finding singletons")
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

    logging.info("Running PBWT matching")

    # There are four params here to mess with:
    #   - query interval
    #   - match group size
    #   - neighbourhood size (L)
    #   - min_matches
    # Keep min_matches low for small sample sizes, can be 2 up to ~1,000 but then > 3
    # Smaller numbers run faster and consume less memory
    MIN_MATCHES = 4
    L = 4
    matcher = Matcher(2 * num_samples, genetic_positions[ac_mask], query_interval, match_group_interval, L, MIN_MATCHES)

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
    logging.info(f"Requested {num_threads} threads, found {actual_num_threads}.")
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
    logging.info(f"Writing to {out}")
    # Save results
    serialize_paths(paths, physical_positions.astype(int), out, start=out_start, end=out_end)
    logging.info(f"Done in (s): {time.time() - start_time}")
    goodbye()

@click.command()
@click.option("--scaffold", required=True, help="Path to vcf containing phased scaffold of common variants")
@click.option("--argn", help="Path to reference ARG in .argn format")
@click.option("--ts", help="Path to reference ARG in .ts format")
@click.option("--unphased", required=True, help="Path to vcf containing the full target dataset (including scaffold variants)")
@click.option("--out", required=True, help="Path to phased output vcf")
def phase(scaffold, argn, ts, unphased, out):
    """
    Use an imputed arg to phase. Other input follows same shape as in SHAPEIT5-rare
    """
    start_time = time.time()
    unphased_vcf = VCF(unphased)
    scaffold_vcf = VCF(scaffold)
    true_vcf = VCF(unphased)
    phased_writer = Writer(out, true_vcf)
    if argn is None and ts is None:
        raise ValueError("Need either --argn or --ts")
    logging.info("Reading ARG...")
    try:
        arg = arg_needle_lib.deserialize_arg(argn)
    except:
        import tskit
        treeseq = tskit.load(ts)
        arg = arg_needle_lib.tskit_to_arg(treeseq)
    arg.populate_children_and_roots()

    i = 0
    logging.info("Phasing...")
    scaffold_empty = False
    v_scaffold = next(scaffold_vcf)
    s_scaffold = 0

    num_hets_found = 0
    # Main phasing routine
    for v in unphased_vcf:
        if not scaffold_empty and v.ID == v_scaffold.ID:
            # If variant exists in scaffold, just copy it
            v.genotypes = v_scaffold.genotypes
            s_scaffold += 1
            try:
                v_scaffold = next(scaffold_vcf)
            except StopIteration:
                scaffold_empty = True
        else:
            # Otherwise, do ARG-based phasing
            G = np.array(v.genotypes)
            g0, g1 = G[:, 0], G[:, 1]
            het_carriers = ((g0 + g1) == 1).nonzero()[0]
            hom_carriers = ((g0 + g1) == 2).nonzero()[0]
            arg_pos = max(v.end - arg.offset, 0)
            arg_pos = min(arg_pos, arg.end - 1)

            for target in het_carriers:
                num_hets_found += 1
                phase_distance_0 = phase_distance(arg, arg_pos, 2 * target, het_carriers, hom_carriers)
                phase_distance_1 = phase_distance(arg, arg_pos, 2 * target + 1, het_carriers, hom_carriers)
                if phase_distance_0 <= phase_distance_1:
                    G[target] = [1, 0, True]
                else:
                    G[target] = [0, 1, True]
            v.genotypes = G

        v.genotypes = v.genotypes
        phased_writer.write_record(v)
    logging.info(f"Done, in {time.time() - start_time} seconds")
    unphased_vcf.close()
    scaffold_vcf.close()
    phased_writer.close()

@click.command()
@click.argument("mode", type=click.Choice(["files", "infer", "convert"]))
@click.option("--threads", required=True, help="Path to an input .threads file.")
@click.option("--argn", default=None, help="Path to an output .argn file.")
@click.option("--tsz", default=None, help="Path to an output .tsz file.")
@click.option("--max_n", default=None, help="How many samples to thread.", type=int)
@click.option("--verify", is_flag=True, show_default=True, default=False, help="Whether to use tskit to verify the ARG.")
def convert(mode, threads, argn, tsz, max_n, verify):
    """
    Convert input .threads file into .threads or .argn file
    """
    start_time = time.time()
    logging.info(f"Starting Threads-{mode} with the following parameters:")
    logging.info(f"  threads: {threads}")
    logging.info(f"  argn:    {argn}")
    logging.info(f"  tsz:     {tsz}")
    logging.info(f"  max_n:   {max_n}")

    if argn is None and tsz is None:
        logging.info("Nothing to do, quitting.")
        sys.exit(0)
    decompressed_threads = decompress_threads(threads)
    try:
        logging.info("Attempting to convert to arg format...")
        arg = threads_to_arg(decompressed_threads, noise=0.0, max_n=max_n, verify=verify)
    except:
        # arg_needle_lib does not allow polytomies
        logging.info(f"Conflicting branches (this is expected), retrying with noise=1e-5...")
        try:
            arg = threads_to_arg(decompressed_threads, noise=1e-5, max_n=max_n, verify=verify)
        except:# tskit.LibraryError:
            logging.info(f"Conflicting branches, retrying with noise=1e-3...")
            arg = threads_to_arg(decompressed_threads, noise=1e-3, max_n=max_n, verify=verify)
    if argn is not None:
        logging.info(f"Writing to {argn}")
        arg_needle_lib.serialize_arg(arg, argn)
    if tsz is not None:
        logging.info(f"Converting to tree sequence and writing to {tsz}")
        tszip.compress(arg_needle_lib.arg_to_tskit(arg), tsz)
    logging.info(f"Done, in {time.time() - start_time} seconds")
    goodbye()

if __name__ == "__main__":
    if len(sys.argv) < 2:
        print("Threads must be called with one of the threads functions: infer, convert, phase.")
        sys.exit()
    else:
        mode = sys.argv[1]
        if mode == "infer":
            infer()
        elif mode == "convert":
            convert()
        elif mode == "phase":
            phase()
        elif mode == "-h" or mode == "--help":
            print("See documentation for each of the Threads functions: infer, convert, phase.")
        else:
            print(f"Unknown mode {mode}")
            sys.exit()
