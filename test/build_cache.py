"""
Pre-generate and cache benchmark data for all sample sizes.

Outputs in threads_cache/:
  sim_{n_dip}.trees          — msprime true tree sequence
  sim_{n_dip}.inferred.trees — tsinfer inferred tree sequence
  sim_{n_dip}.grg_true       — GRG from true ARG
  sim_{n_dip}.grg_inferred   — GRG from inferred ARG
  sim_{n_dip}.genotypes.npz  — genotype matrix + positions
  sim_{n_dip}.threads        — serialized ThreadingInstructions
  sim_{n_dip}.repair         — RePair compressed matrix

Usage:
  python test/build_cache.py [--sizes 50,100,500,...] [--cache-dir threads_cache]
"""

import os, sys, time, argparse, json, glob, resource, platform
import numpy as np
import msprime, tsinfer, pygrgl
import threads_arg, genrepair
from threads_arg.serialization import serialize_instructions, load_instructions

SEED = 42
SEQ_LENGTH = 2e6


def get_rss_bytes():
    """Current resident set size in bytes."""
    ru = resource.getrusage(resource.RUSAGE_SELF)
    if platform.system() == "Darwin":
        return ru.ru_maxrss  # macOS: already bytes
    return ru.ru_maxrss * 1024  # Linux: KB


def measure_memory(fn):
    """Run fn(), return (result, peak_rss_delta_bytes)."""
    import gc; gc.collect()
    before = get_rss_bytes()
    result = fn()
    after = get_rss_bytes()
    return result, max(0, after - before)

def ts_to_biallelic(ts):
    """Extract biallelic sites from tree sequence."""
    positions_bp, positions_cm, genotype_matrix = [], [], []
    for var in ts.variants():
        g = var.genotypes.astype(int).tolist()
        if len(set(g)) < 2 or any(a not in (0, 1) for a in g):
            continue
        positions_bp.append(var.site.position)
        positions_cm.append(var.site.position * 1.3e-8 * 100)
        genotype_matrix.append(g)
    return positions_bp, positions_cm, genotype_matrix


def build_threads_instructions(ts, positions_bp, positions_cm, genotype_matrix):
    """Run full threads pipeline, return ThreadingInstructions."""
    n_haps = ts.num_samples
    matcher = threads_arg.Matcher(n_haps, positions_cm, 0.01, 0.5, 4, 2)
    for g in genotype_matrix:
        matcher.process_site(g)
    target_ids = list(range(1, n_haps))
    match_data = matcher.serializable_matches(target_ids)
    cm_pos = matcher.cm_positions()

    tlm = threads_arg.ThreadsLowMem(
        target_ids, positions_bp, positions_cm, [10000.0], [0.0], 1.4e-8, False)
    tlm.initialize_viterbi(match_data, cm_pos)
    G_np = np.array(genotype_matrix, dtype=np.int32)
    tlm.process_all_sites_viterbi_numpy(G_np)
    tlm.prune()
    tlm.traceback()
    tlm.process_all_sites_hets_numpy(G_np)
    tlm.date_segments()

    all_starts, all_ids, all_heights, all_hetsites = tlm.serialize_paths()
    positions_int = [int(p) for p in positions_bp]
    # serialize_paths returns site indices as starts; convert to physical positions
    all_starts = [[positions_int[s] for s in sample_starts] for sample_starts in all_starts]
    ti = threads_arg.ThreadingInstructions(
        all_starts, all_heights, all_ids, all_hetsites,
        positions_int, positions_int[0], positions_int[-1] + 1)
    return ti



def build_one(n_dip, cache_dir):
    """Generate all cached files for one sample size. Returns timing dict."""
    prefix = os.path.join(cache_dir, f"sim_{n_dip}")
    ts_path = f"{prefix}.trees"
    inf_path = f"{prefix}.inferred.trees"
    grg_true_path = f"{prefix}.grg_true"
    grg_inf_path = f"{prefix}.grg_inferred"
    geno_path = f"{prefix}.genotypes.npz"
    threads_path = f"{prefix}.threads"
    repair_path = f"{prefix}.repair"

    timing = {"n_dip": n_dip}

    print(f"\n{'='*70}")
    print(f"  {n_dip} diploid ({2*n_dip} haploid)")
    print(f"{'='*70}", flush=True)

    # 1. Simulate
    t0 = time.perf_counter()
    ts = msprime.sim_ancestry(
        samples=n_dip, sequence_length=SEQ_LENGTH,
        recombination_rate=1.3e-8, population_size=10000, random_seed=SEED)
    ts = msprime.sim_mutations(ts, rate=1.4e-8, random_seed=SEED + 1)
    ts.dump(ts_path)
    positions_bp, positions_cm, genotype_matrix = ts_to_biallelic(ts)
    n_haps = ts.num_samples
    n_sites = len(genotype_matrix)
    timing["simulate_s"] = time.perf_counter() - t0
    timing["n_hap"] = n_haps
    timing["n_sites"] = n_sites
    print(f"  simulate: {timing['simulate_s']:.1f}s  ({n_haps} hap, {n_sites} sites)", flush=True)

    # 2. Save genotypes
    G_np = np.array(genotype_matrix, dtype=np.int8)
    np.savez_compressed(geno_path,
                        genotypes=G_np,
                        positions_bp=np.array(positions_bp),
                        positions_cm=np.array(positions_cm))
    timing["genotypes_bytes"] = os.path.getsize(geno_path)
    print(f"  genotypes: {timing['genotypes_bytes']/1024:.0f}K", flush=True)

    # 3. GRG from true ARG
    t0 = time.perf_counter()
    grg_true, mem = measure_memory(
        lambda: pygrgl.grg_from_trees(ts_path, binary_mutations=True))
    pygrgl.save_grg(grg_true, grg_true_path)
    timing["grg_true_build_s"] = time.perf_counter() - t0
    timing["grg_true_bytes"] = os.path.getsize(grg_true_path)
    timing["grg_true_mem_bytes"] = mem
    timing["grg_true_mutations"] = grg_true.num_mutations
    print(f"  grg_true: {timing['grg_true_build_s']:.1f}s  mem={mem/1024/1024:.0f}M  "
          f"({grg_true.num_samples} samples, {grg_true.num_mutations} mutations, "
          f"{timing['grg_true_bytes']/1024:.0f}K)", flush=True)

    # 4. tsinfer + GRG from inferred ARG
    t0 = time.perf_counter()
    def _run_tsinfer():
        with tsinfer.SampleData(sequence_length=ts.sequence_length) as sd:
            for var in ts.variants():
                g = var.genotypes
                alleles = var.alleles
                if len(alleles) != 2 or not all(a in (0, 1) for a in g):
                    continue
                if len(set(g)) < 2:
                    continue
                sd.add_site(var.site.position, g, alleles=list(alleles))
        return tsinfer.infer(sd)
    ts_inf, tsinfer_mem = measure_memory(_run_tsinfer)
    ts_inf.dump(inf_path)
    timing["tsinfer_build_s"] = time.perf_counter() - t0
    timing["tsinfer_mem_bytes"] = tsinfer_mem
    print(f"  tsinfer: {timing['tsinfer_build_s']:.1f}s  mem={tsinfer_mem/1024/1024:.0f}M", flush=True)

    t0 = time.perf_counter()
    grg_inf, grg_inf_mem = measure_memory(
        lambda: pygrgl.grg_from_trees(inf_path, binary_mutations=True))
    pygrgl.save_grg(grg_inf, grg_inf_path)
    timing["grg_inf_convert_s"] = time.perf_counter() - t0
    timing["grg_inf_total_s"] = timing["tsinfer_build_s"] + timing["grg_inf_convert_s"]
    timing["grg_inf_bytes"] = os.path.getsize(grg_inf_path)
    timing["grg_inf_mem_bytes"] = tsinfer_mem + grg_inf_mem
    timing["grg_inf_mutations"] = grg_inf.num_mutations
    print(f"  grg_inferred: {timing['grg_inf_convert_s']:.1f}s convert, "
          f"{timing['grg_inf_total_s']:.1f}s total  mem={timing['grg_inf_mem_bytes']/1024/1024:.0f}M  "
          f"({grg_inf.num_samples} samples, {grg_inf.num_mutations} mutations, "
          f"{timing['grg_inf_bytes']/1024:.0f}K)", flush=True)

    # 5. Threads
    t0 = time.perf_counter()
    ti, threads_mem = measure_memory(
        lambda: build_threads_instructions(ts, positions_bp, positions_cm, genotype_matrix))
    timing["threads_build_s"] = time.perf_counter() - t0
    timing["threads_mem_bytes"] = threads_mem
    t0 = time.perf_counter()
    serialize_instructions(ti, threads_path)
    timing["threads_serialize_s"] = time.perf_counter() - t0
    timing["threads_bytes"] = os.path.getsize(threads_path)
    timing["threads_n_samples"] = ti.num_samples
    timing["threads_n_sites"] = ti.num_sites
    print(f"  threads: {timing['threads_build_s']:.1f}s build + "
          f"{timing['threads_serialize_s']:.1f}s serialize  mem={threads_mem/1024/1024:.0f}M  "
          f"({ti.num_samples} hap, {ti.num_sites} sites, "
          f"{timing['threads_bytes']/1024:.0f}K)", flush=True)

    # 7. Prepare tree multiply (measure separately)
    t0 = time.perf_counter()
    _, prep_mem = measure_memory(lambda: ti.prepare_tree_multiply())
    timing["tree_prepare_s"] = time.perf_counter() - t0
    timing["tree_prepare_mem_bytes"] = prep_mem
    print(f"  tree_prepare: {timing['tree_prepare_s']*1000:.1f}ms  mem={prep_mem/1024/1024:.0f}M", flush=True)

    # 6. RePair (saves multiple files with repair_path as basename)
    t0 = time.perf_counter()
    cm, repair_mem = measure_memory(
        lambda: genrepair.CompressedMatrix.from_numpy(
            G_np, repair_path, diploid=False, standardized=False))
    timing["repair_build_s"] = time.perf_counter() - t0
    timing["repair_mem_bytes"] = repair_mem
    repair_files = glob.glob(f"{repair_path}*")
    timing["repair_bytes"] = sum(os.path.getsize(f) for f in repair_files)
    print(f"  repair: {timing['repair_build_s']:.1f}s  mem={repair_mem/1024/1024:.0f}M  "
          f"({timing['repair_bytes']/1024:.0f}K across {len(repair_files)} files)", flush=True)

    all_files = glob.glob(f"{prefix}*")
    total_size = sum(os.path.getsize(f) for f in all_files)
    print(f"  total cached: {total_size/1024/1024:.1f}M ({len(all_files)} files)", flush=True)

    return timing


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--sizes", default="50,100,250,500,1000,2000,5000,10000,20000",
                        help="Comma-separated diploid sample sizes")
    parser.add_argument("--cache-dir", default="threads_cache")
    args = parser.parse_args()

    sizes = [int(s) for s in args.sizes.split(",")]
    os.makedirs(args.cache_dir, exist_ok=True)

    print(f"Building cache for sizes: {sizes}")
    print(f"Cache dir: {args.cache_dir}")

    all_timing = []
    for n_dip in sizes:
        try:
            timing = build_one(n_dip, args.cache_dir)
            all_timing.append(timing)
        except Exception as e:
            print(f"  ERROR at {n_dip} dip: {e}")
            import traceback; traceback.print_exc()

    # Save timing data
    timing_path = os.path.join(args.cache_dir, "build_timing.json")
    with open(timing_path, "w") as f:
        json.dump(all_timing, f, indent=2)
    print(f"\nTiming saved to {timing_path}")

    # Print summary table
    print(f"\n{'='*120}")
    print(f"  CONSTRUCTION TIME & MEMORY SUMMARY")
    print(f"{'='*120}")
    print(f"{'n_dip':>7s} {'hap':>6s} {'sites':>6s} | "
          f"{'threads':>10s} {'tsinfer+grg':>12s} {'RePair':>10s} | "
          f"{'thr_mem':>8s} {'grg_mem':>8s} {'rpr_mem':>8s} | "
          f"{'thr_disk':>9s} {'grg_disk':>9s} {'rpr_disk':>9s}")
    print("-" * 120)
    for t in all_timing:
        def fmt_time(s):
            if s < 60: return f"{s:.1f}s"
            return f"{s/60:.1f}m"
        def fmt_mem(b):
            if b < 1024**2: return f"{b/1024:.0f}K"
            return f"{b/1024/1024:.0f}M"
        def fmt_disk(b):
            if b < 1024**2: return f"{b/1024:.0f}K"
            return f"{b/1024/1024:.1f}M"
        print(f"{t['n_dip']:7d} {t['n_hap']:6d} {t['n_sites']:6d} | "
              f"{fmt_time(t['threads_build_s']):>10s} "
              f"{fmt_time(t['grg_inf_total_s']):>12s} "
              f"{fmt_time(t['repair_build_s']):>10s} | "
              f"{fmt_mem(t['threads_mem_bytes']):>8s} "
              f"{fmt_mem(t['grg_inf_mem_bytes']):>8s} "
              f"{fmt_mem(t['repair_mem_bytes']):>8s} | "
              f"{fmt_disk(t['threads_bytes']):>9s} "
              f"{fmt_disk(t['grg_inf_bytes']):>9s} "
              f"{fmt_disk(t['repair_bytes']):>9s}")


if __name__ == "__main__":
    main()
