"""Simulate ARGs with msprime and infer .threads files at various sample sizes.

Caches results in bench_data/ so bench_multiply.py can load them quickly.
Usage: python simulate_threads.py [n1 n2 ...]
Default sizes: 100 200 500 1000 2000 5000
"""
import os
import sys
import time
import shutil
import tempfile
import numpy as np
import msprime
import pgenlib
from threads_arg.infer import threads_infer

SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
DEMO_FILE = os.path.join(SCRIPT_DIR, "test", "data", "CEU_unscaled.demo")
CACHE_DIR = os.path.join(SCRIPT_DIR, "bench_data")
MUT_RATE = 1.4e-8
RECOMB_RATE = 1.3e-8
SEQ_LENGTH = 5_000_000


def ts_to_pgen(ts, prefix):
    """Write a tskit tree sequence to PGEN/PVAR/PSAM files (biallelic only)."""
    geno = ts.genotype_matrix()  # (num_sites, num_haplotypes)
    # Filter to biallelic sites (drop back-mutations with allele codes > 1)
    biallelic = geno.max(axis=1) <= 1
    geno = geno[biallelic]
    sites = [s for s, keep in zip(ts.sites(), biallelic) if keep]
    num_sites, num_haps = geno.shape
    num_diploid = num_haps // 2

    with open(prefix + ".psam", "w") as f:
        f.write("#IID\tSEX\n")
        for i in range(num_diploid):
            f.write(f"tsk_{i}\tNA\n")

    with open(prefix + ".pvar", "w") as f:
        f.write("#CHROM\tPOS\tID\tREF\tALT\n")
        for i, site in enumerate(sites):
            pos = int(site.position) + 1
            f.write(f"1\t{pos}\tsnp_{i}\tA\tT\n")

    writer = pgenlib.PgenWriter(
        filename=(prefix + ".pgen").encode(),
        sample_ct=num_diploid,
        variant_ct=num_sites,
        hardcall_phase_present=True,
    )
    for s in range(num_sites):
        alleles = np.empty(2 * num_diploid, dtype=np.int32)
        alleles[0::2] = geno[s, 0::2]
        alleles[1::2] = geno[s, 1::2]
        writer.append_alleles(alleles, all_phased=True)
    writer.close()


def simulate_and_infer(n_haps, seed=42):
    """Simulate with msprime, infer with threads, save to cache."""
    out_file = os.path.join(CACHE_DIR, f"n{n_haps}.threads")
    if os.path.exists(out_file):
        print(f"  n={n_haps}: already cached at {out_file}")
        return

    print(f"  n={n_haps}: simulating...", end="", flush=True)
    ts = msprime.sim_ancestry(
        samples=n_haps // 2, ploidy=2, sequence_length=SEQ_LENGTH,
        recombination_rate=RECOMB_RATE, population_size=10_000, random_seed=seed,
    )
    ts = msprime.sim_mutations(ts, rate=MUT_RATE, random_seed=seed)
    print(f" sites={ts.num_sites} trees={ts.num_trees}", end="", flush=True)

    tmpdir = tempfile.mkdtemp(prefix="sim_threads_")
    try:
        prefix = os.path.join(tmpdir, "panel")
        ts_to_pgen(ts, prefix)

        map_file = os.path.join(tmpdir, "gmap.map")
        with open(map_file, "w") as f:
            f.write("pos\tchr\tcM\n")
            f.write(f"0\t1\t0\n")
            f.write(f"{int(SEQ_LENGTH)}\t1\t{SEQ_LENGTH * RECOMB_RATE * 100:.6f}\n")

        print(f" inferring...", end="", flush=True)
        t0 = time.perf_counter()
        threads_infer(
            pgen=prefix + ".pgen",
            map=map_file,
            recombination_rate=RECOMB_RATE,
            demography=DEMO_FILE,
            mutation_rate=MUT_RATE,
            fit_to_data=False,
            allele_ages=None,
            query_interval=0.01,
            match_group_interval=0.5,
            mode="wgs",
            num_threads=1,
            region=None,
            max_sample_batch_size=None,
            save_metadata=False,
            out=out_file,
        )
        elapsed = time.perf_counter() - t0
        print(f" done in {elapsed:.1f}s → {out_file}")
    finally:
        shutil.rmtree(tmpdir, ignore_errors=True)


def main():
    import logging
    logging.basicConfig(level=logging.WARNING)

    os.makedirs(CACHE_DIR, exist_ok=True)

    sizes = [int(x) for x in sys.argv[1:]] if len(sys.argv) > 1 else [100, 200, 500, 1000, 2000, 5000]
    print(f"Simulating + inferring .threads for n = {sizes}")
    print(f"Cache directory: {CACHE_DIR}")

    for n in sizes:
        simulate_and_infer(n)

    print("\nDone. Run bench_multiply.py to benchmark.")


if __name__ == "__main__":
    main()
