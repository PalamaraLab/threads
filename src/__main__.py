import os
import gc
import sys
import time
import h5py
import tszip
import click
import logging
import pgenlib
import subprocess
import arg_needle_lib

import ray
import numpy as np
import xarray as xr
import pandas as pd

from threads import Threads, ThreadsLowMem, Matcher, ViterbiPath
from .utils import decompress_threads, read_map_gz, parse_demography
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
    print(
    """
    `-:-.   ,-;"`-:-.   ,-;"`-:-.   ,-;"`-:-.   ,-;"
       `=`,'=_______________________________=`,'=/
         y== ||Thank you for using Threads|| y==/  
       ,=,-<=‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾=,-<=`.
    ,-'-'   `-=_,-'-'   `-=_,-'-'   `-=_,-'-'   `-=_
    """
    )

def haploidify_samples(samples, out):
    lines = []
    with open(samples, 'r') as shandle:
        for i, l in enumerate(shandle):
            if i <= 1:
                lines.append(l.strip())
            else:
                try:
                    id1, id2, missing = l.strip().split()
                    lines.append(f"{id1}_0 {id2}_0 {missing}")
                    lines.append(f"{id1}_1 {id2}_1 {missing}")
                except:
                    id1, id2, missing, sex = l.strip().split()
                    lines.append(f"{id1}_0 {id2}_0 {missing} {sex}")
                    lines.append(f"{id1}_1 {id2}_1 {missing} {sex}")
    
    with open(out, 'w') as ohandle:
        for l in lines:
            ohandle.write(l + "\n")

def thread_range(haps, bwt, start, end, use_hmm, variant_filter, batch_size, phys_pos):
    threads = []
    range_size = min(end, haps.shape[0])
    # Iterate over genotypes in chunks
    for k in range(int(np.ceil(range_size / batch_size))):
        # Read genotypes into memory
        batch_s = min(end, k * batch_size)
        batch_e = min(end, (k + 1) * batch_size)
        if variant_filter is None:
            haps_chunk = haps[batch_s:batch_e].values.astype(bool)
        else:
            haps_chunk = haps[batch_s:batch_e, variant_filter].values.astype(bool)
        # Iterate over chunk
        for i, g in enumerate(haps_chunk):
            # To save memory, delete hmm when we're done with it
            if use_hmm and batch_s + i >= 100:
                bwt.delete_hmm()
                use_hmm = False

            if i + batch_s < start:
                bwt.insert(g)
            elif i + batch_s == 0:
                bwt.insert(g)
                threads.append(([], [], [], [phys_pos[i] for i, h in enumerate(g) if h]))
            else:
                threads.append(bwt.thread(g))
    return threads

def serialize_paths(paths, positions, out, start=None, end=None):
    # paths = TLM.paths
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
        # bp_ii = path.segment_starts
        # ids = path.sample_ids
        # ages = path.heights
        # het_sites = path.het_sites
        # bps = [positions[bp_idx] for bp_idx in bp_ii]
        if i > 0:
            try:
                assert bps[0] == region_start
                if start is not None:
                    assert bps[0] >= start
                assert len(np.unique(bps)) == len(bps)
                assert bps[-1] <= positions[-1]
            except AssertionError:
                logging.info("funky fragments!")
                logging.info(f"bps[0] {bps[0]}")
                logging.info(f"region_start {region_start}")
                logging.info(f"len(np.unique(bps)) {len(np.unique(bps))}")
                logging.info(f"len(bps) {len(bps)}")
                logging.info(f"bps[-1] {bps[-1]}")
                logging.info(f"positions[-1] {positions[-1]}")
                import pdb
                pdb.set_trace()

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

def serialize(out, threads, samples, positions, arg_range, L=1):
    num_threads = len(samples)
    num_sites = len(positions)

    thread_starts = []
    mut_starts = []
    thread_start = 0
    mut_start = 0
    all_bps, all_ids, all_ages, all_het_sites = [], [], [], []
    for i, thread in enumerate(threads):
        bps, ids, ages, het_sites = thread
        all_bps += bps
        all_ids += ids
        all_ages += ages
        all_het_sites += het_sites
        thread_starts.append(thread_start)
        mut_starts.append(mut_start)
        thread_start += len(bps)
        mut_start += len(het_sites)
    num_stitches = len(all_bps)
    num_mutations = len(all_het_sites)

    assert len(all_bps) == len(all_ids) == len(all_ages)
    assert len(threads) == len(samples)

    f = h5py.File(out, "w")
    f.attrs['datetime_created'] = datetime.now().isoformat()

    compression_opts = 9
    dset_samples = f.create_dataset("samples", (num_threads, 3), dtype=int, compression='gzip',
                                    compression_opts=compression_opts)
    dset_pos = f.create_dataset("positions", (num_sites), dtype=int, compression='gzip',
                                 compression_opts=compression_opts)
    # First L columns are random samples for imputation
    dset_targets = f.create_dataset("thread_targets", (num_stitches, L + 2), dtype=int, compression='gzip',
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

    for i in range(L + 1):
        dset_targets[:, i] = [ids[i] for ids in all_ids]

    dset_targets[:, L + 1] = all_bps

    dset_ages[:] = all_ages

    dset_het_s[:] = all_het_sites

    dset_range[:] = arg_range

    f.close()

# this is only for imputation
def write_vcf(out_vcf_gz, samples, positions, snp_ids, afs, hap1_genotypes, hap2_genotypes, dosages, genotyped_pos):
    chrom = "1"
    # with gzip.open(out_vcf_gz, "wb") as f:
    with open(out_vcf_gz, "w") as f:
        f.write("##fileformat=VCFv4.2\n")
        f.write("##FILTER=<ID=PASS,Description=\"All filters passed\">\n")
        f.write(f"##fileDate={datetime.now().strftime('%m/%d/%Y, %H:%M:%S')}\n")
        f.write("##source=Threads 0.0\n")
        f.write(f"##contig=<ID={chrom}>\n")
        f.write("##FPLOIDY=2\n")
        f.write("##INFO=<ID=IMP,Number=0,Type=Flag,Description=\"Imputed marker\">\n")
        f.write("##INFO=<ID=AF,Number=A,Type=Float,Description=\"Estimated allele frequency\">\n")
        f.write("##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Phased genotypes\">\n")
        f.write("##FORMAT=<ID=DS,Number=A,Type=Float,Description=\"Genotype dosage\">\n")
        f.write(("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t" + "\t".join(samples) + "\n"))
        for (pos, snp_id, af, hap1_var, hap2_var, dosage_var) in zip(positions, snp_ids, afs, hap1_genotypes, hap2_genotypes, dosages):
            imp_str = "" if pos in genotyped_pos else "IMP;"
            gt_strings = [f"{hap_1}|{hap_2}:{dosage:.3f}".rstrip("0").rstrip(".") for hap_1, hap_2, dosage in zip(hap1_var, hap2_var, dosage_var)]
            f.write(("\t".join([chrom, str(pos), snp_id.replace("'", "").replace(":", "_"), "A", "C", ".", "PASS", f"{imp_str}AF={af:.4f}", "GT:DS", "\t".join(gt_strings)]) + "\n"))

def threads_to_arg(thread_dict, noise=0.0, max_n=None):
    threading_instructions = thread_dict['threads']
    pos = thread_dict['positions']
    N = max_n if max_n is not None else np.max(thread_dict['samples'].astype(int)) + 1
    arg_start, arg_end = thread_dict['arg_range']
    # "+ 1" so we can include mutations on the last site... wait are we now adding two many +1's?
    # "+ 2" because.... >:-(
    arg = arg_needle_lib.ARG(0, arg_end - arg_start + 2, reserved_samples=N)
    arg.set_offset(int(arg_start))

    # How this should work:
    for i, t in enumerate(threading_instructions[:N]):
        arg.add_sample(str(i))#, str(i))
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

    # Cycling pass
    # for i, t in enumerate(threading_instructions[N:]):
    #     raise RuntimeError(f"Cycling not currently supported here until arg_needle_lib supports node deletion")
        # arg.lazy_delete(i)
        # arg.add_sample(i, str(i))
        # section_starts, thread_ids, thread_heights = t
        # thread_heights += thread_heights * np.random.normal(0.0, noise, len(thread_heights))
        # arg.thread_sample([s - arg.offset for s in section_starts], thread_ids, thread_heights)

    # tskit is better at checking whether arg makes sense
    logging.info("Verifying ARG...")
    _ = arg_needle_lib.arg_to_tskit(arg)
    return arg

# def parse_vcf(genotypes):
#     raise NotImplementedError

# @click.command()
# @click.argument("mode", type=click.Choice(["files", "infer", "convert"]))
# @click.option("--hap_gz", default=None, help="Path to Oxford-format haps file. Optional for 'files'")
# @click.option("--sample", default=None, help="Path to Oxford-format sample file. Required for 'files' with the hap_gz option")
# @click.option("--vcf", default=None, help="Path to vcf/vcf.gz/bcf file. Optional for 'files'")
# @click.option("--zarr", required=True)
# def files(mode, hap_gz, sample, vcf, zarr):
#     """
#     Prepares funky phased bedfiles for Threads inference as well as the allele count file
#     """
#     # Todo: make this not read everything into memory
#     start_time = time.time()
#     if (hap_gz is None and vcf is None) or zarr is None:
#         print_help_msg(files)
#     if hap_gz is not None and vcf is not None:
#         raise RuntimeError("Please only specify one of hap_gz or vcf")
#     logging.info(f"Starting Threads-{mode} with the following parameters:")
#     if hap_gz is not None:
#         logging.info(f"  hap_gz: {hap_gz}")
#         logging.info(f"  sample: {sample}")
#     else:
#         logging.info(f"  vcf: {vcf}")
#     logging.info(f"  zarr:  {zarr}")

#     if hap_gz is not None:
#         raise NotImplementedError("Not sure if I'll ever support Oxford haps files again. Eugh!")
#     samples = []
#     positions = [v.POS for v in VCF(vcf)]
#     for s in VCF(vcf).samples:
#         samples += [f"{s}_0", f"{s}_1"]
#     N = len(samples)
#     snp_ids = [v.ID for v in VCF(vcf)]
#     M = len(snp_ids)
#     genotypes = np.empty((N, M), dtype=bool)
#     for i, v in enumerate(VCF(vcf)):
#         v_gt = np.array([[g[0], g[1]] for g in v.genotypes]).flatten()
#         genotypes[:, i] = v_gt

#     ds = xr.Dataset(
#         {"genotypes": (("samples", "snp_ids"), genotypes)},
#         coords={
#             "samples": samples,
#             "snp_ids": snp_ids,
#             "pos": np.array(positions),
#             "allele_counts": genotypes.sum(axis=0)
#         }
#     )
#     ds = ds.chunk({"samples": "auto", "snp_ids": -1})
#     ds.to_zarr(zarr, mode="w")
#     end_time = time.time()
#     logging.info(f"Done, in {end_time - start_time} seconds")
#     goodbye()

# @click.command()
# @click.argument("mode", type=click.Choice(["files", "infer", "convert", "impute"]))
# @click.option("--genotypes", required=True, help="zarr thingy")
# @click.option("--map_gz", required=True, help="Path to genotype map")
# @click.option("--adaptive", default=False, is_flag=True, help="If set, ultra-rare sites are gradually trimmed.")
# @click.option("--mutation_rate", required=True, type=float, help="Per-site-per-generation SNP mutation rate.")
# @click.option("--demography", required=True, help="Path to file containing demographic history.")
# @click.option("--modality", required=True, type=click.Choice(['array', 'seq'], case_sensitive=False), help="Inference mode, must be either 'array' or 'seq'.")
# @click.option("--threads", required=True, help="Path to output .threads file.")
# @click.option("--burn_in_left", default=0, help="Left burn-in in base-pairs. Default is 0.")
# @click.option("--burn_in_right", default=0, help="Right burn-in in base-pairs. Default is 0.")
# @click.option("--max_ac", default=5, help="If using adaptive, lowest allele count included in all analyses")
# def infer(mode, genotypes, map_gz, mutation_rate, demography, modality, threads, burn_in_left, burn_in_right, adaptive, max_ac, cycle_n=0):
#     """Infer threading instructions from input data and save to .threads format"""
#     start_time = time.time()
#     logging.info(f"Starting Threads-{mode} with the following parameters:")
#     logging.info(f"  genotypes:      {genotypes}")
#     logging.info(f"  map_gz:         {map_gz}")
#     logging.info(f"  threads:        {threads}")
#     logging.info(f"  demography:     {demography}")
#     logging.info(f"  mutation_rate:  {mutation_rate}")
#     logging.info(f"  cycle_n:        {cycle_n}")
#     logging.info(f"  burn_in_left:   {burn_in_left}")
#     logging.info(f"  burn_in_right:  {burn_in_right}")
#     threads_out = threads

#     if os.path.isdir(genotypes):
#         logging.info(f"{genotypes} is a directory, interpreting as a .zarr archive")
#         z_dataset = xr.open_zarr(genotypes)
#         batch_size = z_dataset.chunks["samples"][0]
#         haps = z_dataset["genotypes"]
#         acounts = z_dataset["allele_counts"].values
#     elif os.path.isfile(genotypes):
#         logging.info(f"{genotypes} is not a .zarr archive, interpreting as a vcf - this will read everything into memory")
#         raise NotImplementedError
#         # This will do the whole conversion from vcf to zarr
#         haps = parse_vcf(genotypes)
#         batch_size = 4e6 // haps.shape[1]
#         acounts = haps.values.sum(axis=0)
#     else:
#         raise ValueError(f"{genotypes} not found")

#     # Keep track of ram
#     max_mem_used_GB = 0

#     cm_pos, phys_pos = read_map_gz(map_gz)
#     logging.info(f"Found a genotype map with {len(phys_pos)} markers")
#     logging.info(f"Found genotype of shape {haps.shape}")

#     phys_pos_0 = phys_pos
#     # arg_range = (phys_pos_0[0] + burn_in_left, phys_pos_0[-1] - burn_in_right + 1)
#     arg_range = (phys_pos_0[0] + burn_in_left, phys_pos_0[-1] - burn_in_right)
#     cm_pos_0 = cm_pos
#     M = len(phys_pos_0)
#     ne_times, ne_sizes = parse_demography(demography)
#     num_samples = haps.shape[0]

#     sparse_sites = modality == "array"
#     use_hmm = modality != "array"
#     threads = []
#     samples = []

#     logging.info("Threading! This could take a while...")
#     if modality == "array" or (modality == "seq" and not adaptive):
#         # batch_size = int(np.ceil(2e6 / M))
#         bwt = Threads(phys_pos, cm_pos, mutation_rate, ne_sizes, ne_times, sparse_sites, use_hmm=use_hmm, burn_in_left=burn_in_left, burn_in_right=burn_in_right)
#         assert (bwt.threading_start, bwt.threading_end) == arg_range
#         for k in range(int(np.ceil(num_samples / batch_size))):
#             s = k * batch_size
#             e = min(num_samples, (k + 1) * batch_size)
#             # TODO: use pandas_plink.Chunk ?
#             haps_chunk = haps[s:e].values.astype(bool)

#             for i, hap in enumerate(haps_chunk):
#                 if (i + s + 1) % 1000 == 0:
#                     logging.info("\n ")

#                 hap = haps_chunk[i]
#                 if (i + s > 0):
#                     thread = bwt.thread(hap)
#                     threads.append(thread)
#                 else:
#                     bwt.insert(hap)
#                     threads.append(([], [], [], [phys_pos[i] for i, h in enumerate(hap) if h]))
#                 samples.append(i + s)
#         if cycle_n > 0:
#             logging.info("Cycling! This shouldn't take too long...")
#             haps_chunk = haps[:min(num_samples, cycle_n)].values.astype(bool)
#             for i, g in enumerate(haps_chunk):
#                 #g = haps[i].values.astype(bool)    
#                 samples.append(i)
#                 bwt.delete_ID(i)
#                 thread = bwt.thread(i, g)
#                 threads.append(thread)
#     else:
#         # counts_path = f"{bfile}.acount"
#         # counts_df = pd.read_table(counts_path)
#         # weird fix because of plink hack
#         # acounts = counts_df['ALT_CTS'].values

#         assert acounts.max() <= num_samples
#         # assert acounts.min() >= num_samples
#         assert len(acounts) == len(phys_pos_0)
#         # acounts -= num_samples

#         bwt = None
#         mumumultiplier = 1.0
#         for ac in range(max_ac + 1):
#             start = int(1 + 1_000 * ac)
#             end = int(1 + 1_000 * (ac + 1))
#             use_hmm = False
#             if ac == 0:
#                 start = 0
#                 use_hmm = True
#                 variant_filter = None
#             else:
#                 if ac == max_ac:
#                     end = num_samples
#                 variant_filter = (ac < acounts) * (acounts < num_samples - ac)
#                 # We still need to get the whole range, though
#                 variant_filter[0] = True
#                 variant_filter[-1] = True  # this one might be unnecessary
#                 phys_pos = phys_pos_0[variant_filter]
#                 cm_pos = cm_pos_0[variant_filter]
#             mumumultiplier = 0.5 * len(phys_pos) / M
#             bwt = Threads(phys_pos, cm_pos, mumumultiplier * mutation_rate, ne_sizes, ne_times, sparse_sites=False, use_hmm=use_hmm, burn_in_left=burn_in_left, burn_in_right=burn_in_right)
#             assert (bwt.threading_start, bwt.threading_end) == arg_range

#             batch_size = int(np.ceil(2e6 / len(phys_pos)))
#             threads += thread_range(haps, bwt, start, end, use_hmm, variant_filter, batch_size, phys_pos)
#             assert len(threads) == bwt.num_samples

#             if end >= num_samples:
#                 break
#         samples = [i for i in range(num_samples)]

#     logging.info(f"Saving results to {threads_out}")
#     serialize(threads_out, threads, samples, phys_pos_0, arg_range)
#     end_time = time.time()
#     logging.info(f"Done, in (s) {end_time - start_time}")
#     goodbye()

# @ray.remote
def partial_viterbi(pgen, mode, num_samples_hap, physical_positions, genetic_positions, demography, mu, sample_batch, s_match_group, match_cm_positions):
    reader = pgenlib.PgenReader(pgen.encode())
    ne_times, ne = parse_demography(demography)

    sparse = None
    if mode == "array":
        sparse = True
    elif mode == "wgs":
        sparse = False
    else:
        raise RuntimeError
    
    TLM = ThreadsLowMem(sample_batch, physical_positions, genetic_positions, ne, ne_times, mu, sparse)
    logging.info("Initializing HMMs")

    TLM.initialize_viterbi(s_match_group, match_cm_positions)
    # TLM.initialize_viterbi(match_subset)
    M = reader.get_variant_ct()
    BATCH_SIZE = int(2e7 // num_samples_hap)
    n_batches = int(np.ceil(M / BATCH_SIZE))

    logging.info("Starting HMM inference")
    a_counter   = 0
    pruned_size = num_samples_hap
    prune_threshold = 10 * num_samples_hap
    prune_count = 0
    last_prune  = 0
    branch_tracker = []
    for b in range(n_batches):
        b_start = b * BATCH_SIZE
        b_end = min(M, (b+1) * BATCH_SIZE)
        g_size = b_end - b_start
        alleles_out = np.empty((g_size, num_samples_hap), dtype=np.int32)
        phased_out = np.empty((g_size, num_samples_hap // 2), dtype=np.uint8)
        reader.read_alleles_and_phasepresent_range(b_start, b_end, alleles_out, phased_out)
        for g in alleles_out:
            # more stable: have hard threshold here and if it's being too restrictive, slightly ease it, x20 seems too spiky
            TLM.process_site_viterbi(g)
            if a_counter % 10 == 0:
                n_branches = TLM.count_branches()
                if n_branches > prune_threshold: #10 * pruned_size:
                    if a_counter - last_prune <= 30:
                        prune_threshold *= 2
                    # logging.info(f"trim at {a_counter}/{M}")
                    TLM.prune()
                    pruned_size = TLM.count_branches()
                    prune_count += 1
                    last_prune = a_counter
                # branch_tracker.append(n_branches)
            a_counter += 1

    logging.info(f"HMMs done (did {prune_count} pruning steps)")
    logging.info("Starting traceback")
    TLM.traceback()

    logging.info("Fetching heterozygous sites")
    for b in range(n_batches):
        b_start = b * BATCH_SIZE
        b_end = min(M, (b+1) * BATCH_SIZE)
        g_size = b_end - b_start
        alleles_out = np.empty((g_size, num_samples_hap), dtype=np.int32)
        phased_out = np.empty((g_size, num_samples_hap // 2), dtype=np.uint8)
        reader.read_alleles_and_phasepresent_range(b_start, b_end, alleles_out, phased_out)
        for g in alleles_out:
            TLM.process_site_hets(g)
    # gc.collect()

    logging.info("Dating segments")
    TLM.date_segments()
    return branch_tracker, TLM.serialize_paths()
    # return TLM.paths

def split_list(list, n):
    """Yield n number of sequential chunks from l."""
    sublists = []
    d, r = divmod(len(list), n)
    for i in range(n):
        si = (d+1)*(i if i < r else r) + d*(0 if i < r else i - r)
        sublists.append(list[si:si+(d+1 if i < r else d)])
    return sublists


# @profile
@click.command()
@click.argument("command", type=click.Choice(["files", "infer", "convert", "impute"]))
@click.option("--pgen", required=True)
@click.option("--map_gz", required=True)
@click.option("--demography", required=True)
@click.option("--mode", required=True, type=click.Choice(['array', 'wgs']))
@click.option("--query_interval", type=float, default=0.01, help="For the preliminary haplotype matching. Given in cM.")
@click.option("--match_group_interval", type=float, default=0.5, help="For the preliminary haplotype matching. Given in cM.")
@click.option("--mutation_rate", required=True, type=float)
@click.option("--num_threads", type=int, default=1)
@click.option("--region", help="Of format 123-345, end-exclusive") 
@click.option("--out")
def infer(command, pgen, map_gz, demography, mutation_rate, query_interval, match_group_interval, mode, num_threads, region, out):
    assert command == "infer"
    start_time = time.time()
    logging.info(f"Starting Threads-infer with the following parameters:")
    logging.info(f"  pgen:                 {pgen}")
    logging.info(f"  map_gz:               {map_gz}")
    logging.info(f"  region:               {region}")
    logging.info(f"  demography:           {demography}")
    logging.info(f"  mutation_rate:        {mutation_rate}")
    logging.info(f"  query_interval:       {query_interval}")
    logging.info(f"  match_group_interval: {match_group_interval}")
    logging.info(f"  num_threads:          {num_threads}")
    logging.info(f"  out:                  {out}")

    genetic_positions, physical_positions = read_map_gz(map_gz)

    out_start, out_end = None, None 
    if region is not None: 
        out_start, out_end = [int(r) for r in region.split("-")]

    reader = pgenlib.PgenReader(pgen.encode())

    num_samples = reader.get_raw_sample_ct()
    num_sites = reader.get_variant_ct()
    logging.info(f"Will build an ARG on {2 * num_samples} haplotypes using {num_sites} sites")

    M = len(physical_positions)
    assert num_sites == M
    BATCH_SIZE = int(4e7 // num_samples)
    n_batches = int(np.ceil(M / BATCH_SIZE))

    logging.info("Finding singletons")
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
        ac_mask[b_start:b_end] = alleles_out.sum(axis=1) > 1

    logging.info("Running PBWT matching")

    # There are four params here to mess with:
    #   - query interval
    #   - match group size
    #   - neighbourhood size (L)
    #   - min_matches
    # NB:
    #   - very important to keep min_matches low for big sample sizes, can be 2 up to 1,000 (?) but then > 3
    #   - might want to consider propagating top k (e.g. 4?) states from prev group instead of just the best
    #   - unclear how much query interval/group size should change with N
    #   - I suspect: raise min_matches to 5 and L to 8
    MIN_MATCHES = 4  # maybe even higher for larger N? or should it scale up?
    L = 4
    matcher = Matcher(2 * num_samples, genetic_positions[ac_mask], query_interval, match_group_interval, L, MIN_MATCHES)  # <- smaller numbers are faster

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
    matcher.propagate_adjacent_matches()

    ### From here we want to parallelise
    actual_num_threads = min(len(os.sched_getaffinity(0)), num_threads)
    logging.info(f"Requested {num_threads} threads, found {actual_num_threads}.")
    paths = []
    branch_trackers = []
    if actual_num_threads > 1:
        sample_batches = split_list(list(range(2 * num_samples)), actual_num_threads)
        s_match_groups = [matcher.serializable_matches(sample_batch) for sample_batch in sample_batches]
        match_cm_positions = matcher.cm_positions()
        del matcher
        del alleles_out
        del phased_out
        gc.collect()
        partial_viterbi_remote = ray.remote(partial_viterbi)
        ray.init()
        results = ray.get([partial_viterbi_remote.remote(
            pgen, 
            mode, 
            2 * num_samples, 
            physical_positions, 
            genetic_positions, 
            demography, 
            mutation_rate, 
            sample_batch, 
            s_match_group, 
            match_cm_positions) for sample_batch, s_match_group in zip(sample_batches, s_match_groups)])
        ray.shutdown()
        for sample_batch, result_tuple in zip(sample_batches, results):
            branch_tracker, result_set = result_tuple
            branch_trackers.append(branch_tracker)
            for sample_id, seg_starts, match_ids, heights, hetsites in zip(sample_batch, *result_set):
                paths.append(ViterbiPath(sample_id, seg_starts, match_ids, heights, hetsites))
    else:
        sample_batch = list(range(2 * num_samples))
        s_match_group = matcher.serializable_matches(sample_batch)
        match_cm_positions = matcher.cm_positions()
        del matcher
        gc.collect()
        branch_tracker, results = partial_viterbi(
            pgen, 
            mode, 
            2 * num_samples, 
            physical_positions, 
            genetic_positions, 
            demography, 
            mutation_rate, 
            sample_batch, 
            s_match_group, 
            match_cm_positions)

        branch_trackers.append(branch_tracker)
        for sample_id, seg_starts, match_ids, heights, hetsites in zip(sample_batch, *results): # type: ignore
            paths.append(ViterbiPath(sample_id, seg_starts, match_ids, heights, hetsites))
    logging.info(f"Writing to {out}, {out}.branches.gz")
    # np.savetxt(f"{out}.branches.gz", np.array(branch_trackers))
    serialize_paths(paths, physical_positions.astype(int), out, start=out_start, end=out_end)
    logging.info(f"Done in (s): {time.time() - start_time}")
    goodbye()

@click.command()
@click.argument("mode", type=click.Choice(["files", "infer", "convert", "impute"]))
@click.option("--panel", required=True, help="Prefix of threads-generated PLINK files for reference panel. Array!!")
@click.option("--panel_wgs", required=True, help="Prefix of threads-generated PLINK files for reference panel. Actual panel!!")
@click.option("--target", required=True, help="Prefix of threads-generated PLINK files for samples to impute.")
@click.option("--map_gz", required=True, help="Path to genotype map")
@click.option("--argn", default=None, help="Path to reference arg")
@click.option("--mutation_rate", required=True, type=float, help="Per-site-per-generation SNP mutation rate.")
@click.option("--demography", required=True, help="Path to file containing demographic history.")
@click.option("--out", required=True, help="Path to output .vcf file.")
@click.option("--start_pos", default=0, help="Where we start imputing")
@click.option("--end_pos", default=0, help="Where we stop imputing")
@click.option("--l", default=128, help="Maximum number of matching states sampled per segment. Default is 128 (too high? too high)")
def impute(mode, panel, panel_wgs, target, map_gz, argn, mutation_rate, demography, out, start_pos, end_pos, l):
    """Thread sequences in bfile_target to sequences in bfile_panel. Assumes data are from genotyping arrays"""
    start_time = time.time()
    L = l
    logging.info(f"Starting Threads-{mode} with the following parameters:")
    logging.info(f"  panel:          {panel}")
    logging.info(f"  panel_wgs:      {panel_wgs}")
    logging.info(f"  target:         {target}")
    logging.info(f"  map_gz:         {map_gz}")
    logging.info(f"  argn:           {argn}")
    logging.info(f"  out:            {out}")
    logging.info(f"  demography:     {demography}")
    logging.info(f"  mutation_rate:  {mutation_rate}")
    logging.info(f"  start_pos:      {start_pos}")
    logging.info(f"  end_pos:        {end_pos}")
    logging.info(f"  L:              {L}")
    vcf_out = out
    RARE_VARIANT_CUTOFF = 0.01

    if os.path.isdir(panel):
        logging.info(f"{panel} is a directory, interpreting as a .zarr archive")
        panel_dataset = xr.open_zarr(panel)
        batch_size = panel_dataset.chunks["samples"][0]
        haps_panel = panel_dataset["genotypes"]
    elif os.path.isfile(panel):
        logging.info(f"{panel} is not a .zarr archive, interpreting as a vcf - this will read everything into memory")
        raise NotImplementedError
        # This will do the whole conversion from vcf to zarr
        haps_panel = parse_vcf(panel)
        batch_size = 4e6 // haps_panel.shape[1]
    else:
        raise ValueError(f"{panel} not found")

    if os.path.isdir(target):
        logging.info(f"{target} is a directory, interpreting as a .zarr archive")
        target_dataset = xr.open_zarr(target)
        haps_target = target_dataset["genotypes"]
    elif os.path.isfile(target):
        logging.info(f"{target} is not a .zarr archive, interpreting as a vcf - this will read everything into memory")
        raise NotImplementedError
        # This will do the whole conversion from vcf to zarr
        haps_target = parse_vcf(target)
    else:
        raise ValueError(f"{target} not found")

    # Keep track of ram
    max_mem_used_GB = 0

    cm_pos, phys_pos = read_map_gz(map_gz)
    genotyped_pos = [int(p) for p in phys_pos]
    logging.info(f"Found a genotype map with {len(phys_pos)} markers")

    # haps_panel = read_plink1_bin(f"{bfile_panel}.bed", f"{bfile_panel}.bim", f"{bfile_panel}.fam")
    # haps_target = read_plink1_bin(f"{bfile_target}.bed", f"{bfile_target}.bim", f"{bfile_target}.fam")
    if haps_panel.shape[1] != haps_target.shape[1]:
        raise ValueError(f"Panel has shape {haps_panel.shape} and target has shape {haps_target.shape}")
    elif haps_panel.shape[0] == 0:
        raise ValueError("No sequences found in panel.")
    elif haps_target.shape[0] == 0:
        raise ValueError("No sequences found in target.")

    logging.info(f"Found panel of shape {haps_panel.shape}")
    logging.info(f"Found target of shape {haps_target.shape}")

    M = len(phys_pos)
    ne_times, ne_sizes = parse_demography(demography)
    num_samples_panel = haps_panel.shape[0]
    num_samples_target = haps_target.shape[0]

    sparse_sites = True
    use_hmm = False

    bwt = Threads(phys_pos,
                  cm_pos,
                  mutation_rate,
                  ne_sizes,
                  ne_times,
                  sparse_sites,
                  use_hmm=use_hmm)

    logging.info("Building panel...")
    for k in range(int(np.ceil(num_samples_panel / batch_size))):
        s = k * batch_size
        e = min(num_samples_panel, (k + 1) * batch_size)
        haps_chunk = haps_panel[s:e].values.astype(bool)
        for i, hap in enumerate(haps_chunk):
            hap = haps_chunk[i]
            bwt.insert(hap)

    logging.info("Done building panel, imputing...")
    imputation_threads = []
    for k in range(int(np.ceil(num_samples_target / batch_size))):
        s = k * batch_size
        e = min(num_samples_target, (k + 1) * batch_size)
        haps_chunk = haps_target[s:e].values.astype(bool)

        for i, hap in enumerate(haps_chunk):
            hap = haps_chunk[i]
            imputation_thread = bwt.impute(hap, L)
            imputation_threads.append(imputation_thread)
            # print(len(imputation_thread))
    assert len(imputation_threads) == num_samples_target
    logging.info(f"Done HMM stuff, now doing all the dirty copying work")

    logging.info(f"{panel_wgs} is a directory, interpreting as a .zarr archive")
    panel_dataset = xr.open_zarr(panel_wgs)
    batch_size = panel_dataset.chunks["samples"][0]
    positions = panel_dataset["pos"].values
    snp_ids = panel_dataset["snp_ids"]
    # incl_idx = (start_pos <= positions) * (positions <= end_pos)
    haps = panel_dataset["genotypes"]
    acount = panel_dataset["allele_counts"].values
    M = haps.shape[1]
    mafs = acount / haps.shape[0]
    panel_ids = panel_dataset["samples"].values#haps.iid.values

    mut_lookups = 0
    dosage_lookups = 0

    # Get actual sample names, this is an ugly workaround until we actually save sample names with threads (lazy)
    N_dip = haps.shape[0] // 2
    diploid_samples = []
    for k in range(N_dip, N_dip + num_samples_target // 2):
        diploid_samples.append(f"tsk_{k}")

    arg = None
    if argn is not None:
        arg = arg_needle_lib.deserialize_arg(argn)
        arg.populate_mutations_on_edges()
        # an overabundance of caution here
        arg_sites = (np.array(arg.get_sites()) + arg.offset).astype(int)
        assert len(arg_sites) == M == len(mafs)
        assert (arg_sites == positions).all()

    haps = haps.values
    imputed_genotypes = np.empty((num_samples_target, M), dtype=float)
    imputed_genotypes[:] = np.nan
    # eugh, clean this up, it's gotten just as ugly as the other thing
    for target_id, imputation_thread in enumerate(imputation_threads):
        num_segs = len(imputation_thread)
        num_imputed = 0
        for i, imp_segment in enumerate(imputation_thread):
            seg_end = end_pos + 1 if i == num_segs - 1 else imputation_thread[i + 1].seg_start
            if seg_end <= start_pos:
                continue
            if i == 0 or imp_segment.seg_start < start_pos:
                seg_start = start_pos
            else:
                seg_start = imp_segment.seg_start
            if seg_start > end_pos:
                break
            pos_indices = (seg_start <= positions) * (positions < seg_end)
            seg_positions = positions[pos_indices]
            M_batch = len(seg_positions)
            if M_batch == 0:
                continue
            seg_start_index = np.nonzero(pos_indices)[0][0]#max(, num_imputed)
            if num_imputed > seg_start_index:
                raise RuntimeError
                # this whole things is pretty nuts and shouldn't even be necessary, should take care of this c++-side
                # will only happen with funky impute-overlaps (currently switched off)
                offset = (num_imputed - seg_start_index)
                M_batch -= offset
                seg_positions = seg_positions[offset:]
                seg_start_index = num_imputed
                pos_indices[:seg_start_index] = False
            if M_batch <= 0:
                continue
            # if i == num_segs - 1 or imp_segment.ids != imputation_thread[i + 1].ids:
            # we're between overlap regions and we just average
            seg_ids = np.array(imp_segment.ids)  # [np.array(imp_segment.weights) != 0]
            # current chunking is hilariously bad for this
            seg_haps = haps[seg_ids][:, pos_indices]#.values
            assert np.isnan(imputed_genotypes[target_id, seg_start_index:seg_start_index + M_batch]).all()
            imputed_genotypes[target_id, seg_start_index:seg_start_index + M_batch] = seg_haps.mean(axis=0)
            # else:
            #     next_segment = imputation_thread[i + 1]
            #     # now we're going to interpolate
            #     start_weights = np.array(imp_segment.weights)
            #     end_weights = np.array(next_segment.weights)
            #     seg_haps = haps[imp_segment.ids][:, pos_indices]#.values
            #     pos_interpolation = (seg_positions - seg_start) / (seg_end - seg_start)
            #     assert 0 <= pos_interpolation.min() and pos_interpolation.max() <= 1
            #     seg_weights = np.matmul(start_weights[:, None], 1 - pos_interpolation[None, :]) + np.matmul(end_weights[:, None], pos_interpolation[None, :])
            #     assert np.abs(seg_weights.sum(axis=0) - 1).max() < 0.001
            #     assert seg_weights.shape == seg_haps.shape
            #     assert np.isnan(imputed_genotypes[target_id, seg_start_index:seg_start_index + M_batch]).all()
            #     imputed_genotypes[target_id, seg_start_index:seg_start_index + M_batch] = (seg_weights * seg_haps).sum(axis=0)
            num_imputed += M_batch
            if arg is not None:
                # This is the imputation bit
                mut_indices = np.nonzero(haps[imp_segment.ids, seg_start_index:seg_start_index + M_batch].sum(axis=0))[0]
                mut_mafs = mafs[seg_start_index + mut_indices]
                for mi, mut_maf in zip(mut_indices, mut_mafs):
                    if mut_maf >= RARE_VARIANT_CUTOFF:
                        continue
                    dosages = []
                    # for ti_c, tbi_c in zip(ti, tbi):
                    for p_id, height in zip(imp_segment.ids, imp_segment.ages):  # ti_c, tbi_c in np.unique([(a, b) for (a, b) in zip(ti, tbi)], axis=0):
                        if haps[p_id, seg_start_index + mi] == 1:
                            # arguments are: arg, sample, site, height
                            # print(ti_c, mi+segment_offset, mut_maf)
                            try:
                                _ = arg_needle_lib.most_recent_mutation(arg, p_id, mi + seg_start_index)
                            except RuntimeError:
                                site_muts = [m for m in arg.mutations() if m.site_id == mi + seg_start_index]
                                import pdb
                                pdb.set_trace()
                            dosages.append(arg_needle_lib.threading_dosage(arg, p_id, mi + seg_start_index, height, Ne=40_000))
                            dosage_lookups += 1
                        else:
                            dosages.append(0.)
                    imputed_genotypes[target_id, seg_start_index + mi] = np.mean(dosages)

                # raise NotImplementedError
        try:
            assert num_imputed == M
        except AssertionError:
            import pdb
            pdb.set_trace()
            # TODO: the arg bit
            # imputed_genotypes += M_batch - seg_start_index
        # assert imputed_genotypes == M
    assert np.isnan(imputed_genotypes).sum() == 0
    # merge haps to dips
    diploid_dosages = imputed_genotypes[::2] + imputed_genotypes[1::2]
    hap1_genotypes = np.rint(imputed_genotypes[::2]).astype(int)
    hap2_genotypes = np.rint(imputed_genotypes[1::2]).astype(int)

    # diploid_samples = [sid.replace("_0", "") for sid in target_sample_names[::2]]

    logging.info(f"Writing imputed haplotypes to {out}")
    write_vcf(out, diploid_samples, list(positions), list(snp_ids.values), mafs, hap1_genotypes.transpose(), hap2_genotypes.transpose(), diploid_dosages.transpose(), genotyped_pos)
    end_time = time.time()
    logging.info(f"Done, in (s) {end_time - start_time}")
    goodbye()

@click.command()
@click.argument("mode", type=click.Choice(["phase"]))
@click.option("--panel", required=True)
@click.option("--target", required=True)
@click.option("--out", required=True)
@click.option("--demography", required=True)
@click.option("--mutation_rate", default=1.65e-8, type=float, help="Per-site-per-generation SNP mutation rate.")
@click.option("--map_gz", required=True)
@click.option("--argn", default=None)
def phase(mode, panel, target, out, demography, map_gz, argn, mutation_rate):
    if os.path.isdir(panel):
        logging.info(f"{panel} is a directory, interpreting as a .zarr archive")
        panel_dataset = xr.open_zarr(panel)
        batch_size = panel_dataset.chunks["samples"][0]
        haps_panel = panel_dataset["genotypes"]
    elif os.path.isfile(panel):
        logging.info(f"{panel} is not a .zarr archive, interpreting as a vcf - this will read everything into memory")
        raise NotImplementedError
    else:
        raise ValueError(f"{panel} not found")

    if os.path.isdir(target):
        logging.info(f"{target} is a directory, interpreting as a .zarr archive")
        target_dataset = xr.open_zarr(target)
        dips_target = target_dataset["genotypes"]
    elif os.path.isfile(target):
        dips_target = np.array([[g[0] + g[1] for g in v.genotypes] for v in VCF(target)]).transpose()
    else:
        raise ValueError(f"{target} not found")
    
    cm_pos, phys_pos = read_map_gz(map_gz)
    genotyped_pos = [int(p) for p in phys_pos]
    logging.info(f"Found a genotype map with {len(phys_pos)} markers")

    # haps_panel = read_plink1_bin(f"{bfile_panel}.bed", f"{bfile_panel}.bim", f"{bfile_panel}.fam")
    # haps_target = read_plink1_bin(f"{bfile_target}.bed", f"{bfile_target}.bim", f"{bfile_target}.fam")
    if haps_panel.shape[1] != dips_target.shape[1]:
        raise ValueError(f"Panel has shape {haps_panel.shape} and target has shape {dips_target.shape}")
    elif haps_panel.shape[0] == 0:
        raise ValueError("No sequences found in panel.")
    elif dips_target.shape[0] == 0:
        raise ValueError("No sequences found in target.")

    logging.info(f"Found panel of shape {haps_panel.shape}")
    logging.info(f"Found target of shape {dips_target.shape}")

    M = len(phys_pos)
    ne_times, ne_sizes = parse_demography(demography)
    num_samples_panel = haps_panel.shape[0]
    num_samples_target = dips_target.shape[0]

    sparse_sites = False
    use_hmm = False

    bwt = Threads(phys_pos,
                  cm_pos,
                  mutation_rate,
                  ne_sizes,
                  ne_times,
                  sparse_sites,
                  use_hmm=use_hmm)

    logging.info("Building panel...")
    for k in range(int(np.ceil(num_samples_panel / batch_size))):
        s = k * batch_size
        e = min(num_samples_panel, (k + 1) * batch_size)
        haps_chunk = haps_panel[s:e].values.astype(bool)
        for i, hap in enumerate(haps_chunk):
            hap = haps_chunk[i]
            bwt.insert(hap)
        #     if bwt.num_samples > 10:
        #         break
        # if bwt.num_samples > 10:
        #     break

    logging.info("Phasing...")
    haps_a = np.empty((num_samples_target, M))
    haps_b = np.empty((num_samples_target, M))
    for j, dip in enumerate(dips_target):
        hap_a, hap_b = bwt.phase(dip)
        haps_a[j] = hap_a
        haps_b[j] = hap_b

    logging.info("Writing...")
    # import pdb
    # pdb.set_trace()
    vcf = VCF(target)
    w = Writer(out, vcf)
    for v, hap_a, hap_b in zip(vcf, haps_a.transpose(), haps_b.transpose()):
        for i in range(len(v.genotypes)):
            v.genotypes[i] = [hap_a[i], hap_b[i], True]
            # v.genotypes[i][1] = False
        # need to reassign for some reason(?)
        v.genotypes = v.genotypes
        w.write_record(v)
    w.close()
    vcf.close()

@click.command()
@click.argument("mode", type=click.Choice(["files", "infer", "convert"]))
@click.option("--threads", required=True, help="Path to an input .threads file.")
@click.option("--argn", default=None, help="Path to an output .argn file.")
@click.option("--tsz", default=None, help="Path to an output .tsz file.")
@click.option("--max_n", default=None, help="Path to an output .tsz file.", type=int)
def convert(mode, threads, argn, tsz, max_n):
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
        arg = threads_to_arg(decompressed_threads, noise=0.0, max_n=max_n)
    except:# tskit.LibraryError:
        logging.info(f"Conflicting branches, retrying with noise=1e-5...")
        try:
            arg = threads_to_arg(decompressed_threads, noise=1e-5, max_n=max_n)
        except:# tskit.LibraryError:
            logging.info(f"Conflicting branches, retrying with noise=1e-3...")
            arg = threads_to_arg(decompressed_threads, noise=1e-3, max_n=max_n)
    if argn is not None:
        logging.info(f"Writing to {argn}")
        arg_needle_lib.serialize_arg(arg, argn)
    if tsz is not None:
        logging.info(f"Writing to {tsz}")
        tszip.compress(arg_needle_lib.arg_to_tskit(arg), tsz)
    logging.info(f"Done, in {time.time() - start_time} seconds")
    goodbye()

if __name__ == "__main__":
    if len(sys.argv) < 2:
        print("Threads must be called with one of the threads functions: files, infer, convert, impute.")
        sys.exit()
    else:
        mode = sys.argv[1]
        if mode == "files":
            files()
        elif mode == "infer":
            infer()
        elif mode == "convert":
            convert()
        elif mode == "impute":
            impute()
        elif mode == "phase":
            phase()
        elif mode == "-h" or mode == "--help":
            print("See documentation for each of the Threads functions: files, infer, convert, impute.")
        else:
            print(f"Unknown mode {mode}")
            sys.exit()

