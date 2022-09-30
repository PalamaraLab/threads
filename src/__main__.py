import sys
import time
import h5py
import tszip
import click
import logging
import subprocess
import arg_needle_lib

import numpy as np
import pandas as pd

from threads import DPBWT
from datetime import datetime
from pandas_plink import read_plink1_bin

logging.basicConfig(
    format='%(asctime)s %(levelname)-8s %(message)s',
    level=logging.INFO,
    datefmt='%Y-%m-%d %H:%M:%S')

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

def read_map_gz(map_gz):
    """
        Reading in haps and maps file for Li-Stephens
    """
    if (map_gz[:-3] == ".gz") :
        maps = pd.read_table(map_gz, header=None, compression='gzip')
    else:
        maps = pd.read_table(map_gz, header=None)
    cm_pos = maps[2].values.astype(np.float64)
    phys_pos = maps[3].values.astype(np.float64)
    for i in range(1, len(cm_pos)):
        if cm_pos[i] <= cm_pos[i-1]:
            cm_pos[i] = cm_pos[i-1] + 1e-5
    return cm_pos, phys_pos

def parse_demography(demography):
    d = pd.read_table(demography, delim_whitespace=True, header=None)
    return list(d[0]), list(d[1])

def thread_range(haps, bwt, start, end, use_hmm, variant_filter, batch_size):
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
                threads.append(([], [], []))
            else:
                threads.append(bwt.thread(g))
    return threads

def serialize(out, threads, samples, positions):
    num_threads = len(samples)
    num_sites = len(positions)

    thread_starts = []
    thread_start = 0
    all_bps, all_ids, all_ages = [], [], []
    for i, thread in enumerate(threads):
        # try:
        if i == 0:
            bps, ids, ages = [], [], []
        else:
            bps, ids, ages = thread
        all_bps += bps
        all_ids += ids
        all_ages += ages
        thread_starts.append(thread_start)
        thread_start += len(bps)
    num_stitches = len(all_bps)

    assert len(all_bps) == len(all_ids) == len(all_ages)
    assert len(threads) == len(samples)

    f = h5py.File(out, "w")
    f.attrs['datetime_created'] = datetime.now().isoformat()

    compression_opts = 9
    dset_samples = f.create_dataset("samples", (num_threads, 2), dtype=int, compression='gzip',
                                    compression_opts=compression_opts)
    dset_pos = f.create_dataset("positions", (num_sites), dtype=np.double, compression='gzip',
                                 compression_opts=compression_opts)
    dset_targets = f.create_dataset("thread_targets", (num_stitches, 2), dtype=int, compression='gzip',
                                    compression_opts=compression_opts)
    dset_ages = f.create_dataset("thread_ages", (num_stitches), dtype=np.double, compression='gzip',
                                 compression_opts=compression_opts)
    dset_samples[:, 0] = samples
    dset_samples[:, 1] = thread_starts
    dset_pos[:] = positions

    dset_targets[:, 0] = all_ids
    dset_targets[:, 1] = all_bps

    dset_ages[:] = all_ages

    f.close()

def threads_to_arg(thread_dict, noise=0.0):
    threading_instructions = thread_dict['threads']
    pos = thread_dict['positions']
    N = np.max(thread_dict['samples'].astype(int)) + 1
    arg = arg_needle_lib.ARG(0, pos[-1] - pos[0], reserved_samples=N)
    arg.set_offset(int(pos[0]))

    # How this should work:
    for i, t in enumerate(threading_instructions[:N]):
        arg.add_sample(str(i))#, str(i))
        if i > 0:
            section_starts, thread_ids, thread_heights = t
            thread_heights += thread_heights * np.random.normal(0.0, noise, len(thread_heights))
            # this is weird, should fix
            arg_starts = [s - arg.offset for s in section_starts]
            if arg_starts[-1] >= arg.end:
                arg.thread_sample([s - arg.offset for s in section_starts[:-1]], thread_ids[:-1], thread_heights[:-1])
            else:
                arg.thread_sample([s - arg.offset for s in section_starts], thread_ids, thread_heights)

    # Cycling pass
    for i, t in enumerate(threading_instructions[N:]):
        raise RuntimeError(f"Cycling not currently supported here until we get node deletion into arg_needle_lib")
        arg.lazy_delete(i)
        arg.add_sample(i, str(i))
        section_starts, thread_ids, thread_heights = t
        thread_heights += thread_heights * np.random.normal(0.0, noise, len(thread_heights))
        arg.thread_sample([s - arg.offset for s in section_starts], thread_ids, thread_heights)
    # tskit is better at checking whether arg makes sense
    _ = arg_needle_lib.arg_to_tskit(arg)
    return arg

def decompress_threads(threads):
    f = h5py.File(threads, "r")

    samples, thread_starts = f["samples"][:, 0], f["samples"][:, 1] #dset_samples[:, 0], dset_ = f['flags'][...]
    positions = f['positions'][...]
    flat_ids, flat_bps = f['thread_targets'][:, 0], f['thread_targets'][:, 1]
    flat_ages = f['thread_ages'][...]

    threading_instructions = []
    for i, start in enumerate(thread_starts):
        if i == len(thread_starts) - 1:
            ids = flat_ids[start:]
            bps = flat_bps[start:]
            ages = flat_ages[start:]
        else:
            ids = flat_ids[start:thread_starts[i + 1]]
            bps = flat_bps[start:thread_starts[i + 1]]
            ages = flat_ages[start:thread_starts[i + 1]]
        threading_instructions.append((bps, ids, ages))
    return {
        "threads": threading_instructions,
        "samples": samples,
        "positions": positions
    }

def files(hap_gz, sample, bfile):
    """
    Prepares funky phased bedfiles for Threads inference as well as the allele count file
    """
    start_time = time.time()
    if hap_gz is None:
        raise ValueError("Need to specify a hap_gz")
    if sample is None:
        raise ValueError("Need to specify a sample")
    if bfile is None:
        raise ValueError("Need to specify a bfile")
    logging.info("Starting Threads-file preparation with the following parameters:")
    logging.info(f"  hap_gz: {bfile}")
    logging.info(f"  sample: {sample}")
    logging.info(f"  bfile:  {bfile}")

    new_sample = f"{bfile}.new.sample"
    new_hap_gz = f"{bfile}.new.hap.gz"

    logging.info(f"Making tmp file {new_sample}")
    haploidify_samples(sample, new_sample)
    # This is ripe for injection
    logging.info(f"Making tmp file {new_hap_gz}")
    # This adds an empty column next to all the haps to trick plink1 into reading haploids
    awk_cmd = '{printf "1 "; for (i=2; i<6; i++) printf $i " "; for (i=6; i<NF; i++) printf $i " 0 "; print $NF " 0"}'
    subprocess.run(f"zcat {hap_gz} | awk '{awk_cmd}' | gzip > {new_hap_gz}", shell=True, check=True)
    subprocess.run(["plink2", "--make-bed", "--haps", new_hap_gz, "--sample", new_sample, "--freq", "counts", "--out", bfile], check=True)
    logging.info(f"Removing tmp file {new_sample}")
    subprocess.run(f"rm {new_sample}", shell=True, check=True)
    logging.info(f"Removing tmp file {new_hap_gz}")
    subprocess.run(f"rm {new_hap_gz}", shell=True, check=True)
    end_time = time.time()
    logging.info(f"Done, in {end_time - start_time} seconds")

def infer(bfile, map_gz, mutation_rate, demography, modality, threads, adaptive=False, max_ac=5, cycle_n=0):
    start_time = time.time()
    logging.info("Starting Threads-inference with the following parameters:")
    logging.info(f"  bfile:          {bfile}")
    logging.info(f"  map_gz:         {map_gz}")
    logging.info(f"  threads:        {threads}")
    logging.info(f"  demography:     {demography}")
    logging.info(f"  mutation_rate:  {mutation_rate}")
    logging.info(f"  cycle_n:        {cycle_n}")
    logging.info(f"  modality:       {modality}")
    if bfile is None:
        raise RuntimeError("Missing required argument --bfile")
    if map_gz is None:
        raise RuntimeError("Missing required argument --map_gz")
    if threads is None:
        raise RuntimeError("Missing required argument --threads")
    if bfile is None:
        raise RuntimeError("Missing required argument --bfile")
    if demography is None:
        raise RuntimeError("Missing required argument --demography")
    if mutation_rate is None:
        raise RuntimeError("Missing required argument --mutation_rate")
    threads_out = threads

    cm_pos, phys_pos = read_map_gz(map_gz)
    logging.info(f"Found a genotype map with {len(phys_pos)} markers")
    haps = read_plink1_bin(f"{bfile}.bed", f"{bfile}.bim", f"{bfile}.fam")
    logging.info(f"Found genotype of shape {haps.shape}")

    phys_pos_0 = phys_pos
    cm_pos_0 = cm_pos
    M = len(phys_pos_0)
    ne_times, ne_sizes = parse_demography(demography)
    num_samples = haps.shape[0]

    sparse_sites = modality == "snp"
    use_hmm = modality == "snp"
    threads = []
    samples = []

    logging.info("Threading! This could take a while...")
    if modality == "snp" or (modality == "seq" and not adaptive):
        batch_size = int(np.ceil(2e6 / M))
        bwt = DPBWT(phys_pos, cm_pos, mutation_rate, ne_sizes, ne_times, sparse_sites, use_hmm=use_hmm)
        for k in range(int(np.ceil(num_samples / batch_size))):
            s = k * batch_size
            e = min(num_samples, (k + 1) * batch_size)
            haps_chunk = haps[s:e].values.astype(bool)

            for i, hap in enumerate(haps_chunk):
                if (i + s + 1) % 1000 == 0:
                    logging.info("\n ")

                hap = haps_chunk[i]
                if (i + s > 0):
                    thread = bwt.thread(hap)
                    threads.append(thread)
                else:
                    bwt.insert(hap)
                    threads.append(([], [], []))
                samples.append(i + s)
        if cycle_n > 0:
            logging.info("Cycling! This shouldn't take too long...")
            haps_chunk = haps[:min(num_samples, cycle_n)].values.astype(bool)
            for i, g in enumerate(haps_chunk):
                #g = haps[i].values.astype(bool)    
                samples.append(i)
                bwt.delete_ID(i)
                thread = bwt.thread(i, g)
                threads.append(thread)
    else:
        counts_path = f"{bfile}.acount"
        counts_df = pd.read_table(counts_path)
        # weird fix because of plink hack
        acounts = counts_df['ALT_CTS'].values
        assert acounts.min() >= num_samples
        assert len(acounts) == len(phys_pos_0)
        acounts -= num_samples

        bwt = None
        mumumultiplier = 1.0
        for ac in range(max_ac + 1):
            start = int(1 + 1_000 * ac)
            end = int(1 + 1_000 * (ac + 1))
            # use_hmm = False
            if ac == 0:
                start = 0
                use_hmm = True
                variant_filter = None
            else:
                use_hmm = False
                if ac == max_ac:
                    end = num_samples
                variant_filter = (ac < acounts) & (acounts < num_samples - ac)
                # We still need to get the whole range, though
                variant_filter[0] = True
                variant_filter[-1] = True  # this one might be unnecessary
                phys_pos = phys_pos_0[variant_filter]
                cm_pos = cm_pos_0[variant_filter]
            mumumultiplier = 0.5 * len(phys_pos) / M
            bwt = DPBWT(phys_pos, cm_pos, mumumultiplier * mutation_rate, ne_sizes, ne_times, sparse_sites=False, use_hmm=use_hmm)

            batch_size = int(np.ceil(2e6 / len(phys_pos)))
            threads += thread_range(haps, bwt, start, end, use_hmm, variant_filter, batch_size)
            assert len(threads) == bwt.num_samples

            if end >= num_samples:
                break
        samples = [i for i in range(num_samples)]

    logging.info(f"Saving results to {threads_out}")
    serialize(threads_out, threads, samples, phys_pos_0)
    end_time = time.time()
    logging.info(f"Done, in (s): {end_time - start_time}")

def convert(threads, argn, tsz):
    start_time = time.time()
    logging.info("Starting Threads-inference with the following parameters:")
    logging.info(f"  threads: {threads}")
    logging.info(f"  argn:    {argn}")
    logging.info(f"  tsz:     {tsz}")
    if threads is None:
        logging.info(f"Missing required argument --threads")
    if argn is None and tsz is None:
        logging.info(f"Nothing to do, quitting.")
        sys.exit(0)
    decompressed_threads = decompress_threads(threads)
    try:
        arg = threads_to_arg(decompressed_threads, noise=0.0)
    except:# tskit.LibraryError:
        try:
            arg = threads_to_arg(decompressed_threads, noise=1e-5)
        except:# tskit.LibraryError:
            arg = threads_to_arg(decompressed_threads, noise=1e-3)
    if argn is not None:
        logging.info(f"Writing to {argn}")
        arg_needle_lib.serialize_arg(arg, argn)
    if tsz is not None:
        logging.info(f"Writing to {tsz}")
        tszip.compress(arg_needle_lib.arg_to_tskit(arg), tsz)
    logging.info(f"Done, in {time.time() - start_time} seconds")

@click.command()
@click.argument("mode", type=click.Choice(["files", "infer", "convert"]))
@click.option("--hap_gz", default=None, help="required for 'files'")
@click.option("--sample", default=None, help="required for 'files'")
@click.option("--bfile", default=None, help="required for 'files' and 'infer")
@click.option("--map_gz", default=None, help="required for 'infer'")
@click.option("--adaptive", default=False, is_flag=True, help="optional for 'infer'")
@click.option("--mutation_rate", default=None, type=float, help="required for 'infer'")
@click.option("--demography", default=None, help="required for 'infer'")
@click.option("--modality", default=None, help="required for 'infer'")
@click.option("--threads", default=None, help="required for 'infer' and 'convert'")
@click.option("--argn", default=None, help="required for 'convert'")
@click.option("--tsz", default=None, help="required for 'convert'")
def main(mode, hap_gz, sample, bfile, map_gz, adaptive, mutation_rate, demography, modality, threads, argn, tsz):
    if mode == "files":
        files(hap_gz, sample, bfile)
    elif mode == "infer":
        infer(bfile, map_gz, mutation_rate, demography, modality, threads, adaptive=adaptive)
    elif mode == "convert":
        convert(threads, argn, tsz)
    else:
        raise ValueError(f"Unknown mode {mode}")

if __name__ == "__main__":
    main()
