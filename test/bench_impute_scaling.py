#!/usr/bin/env python3
"""
Scaling microbenchmark for imputation hot paths.

Isolates fwbw, sparsify, per-variant dosage computation, posterior cache
rebuild, and VCF write formatting from VCF I/O.
Uses msprime to generate synthetic data at configurable scale.

Usage:
    python test/bench_impute_scaling.py
    python test/bench_impute_scaling.py --preset medium
    python test/bench_impute_scaling.py --preset large
    python test/bench_impute_scaling.py --n-panel 2000 --n-target 100
"""
import argparse
import time
import sys

import numpy as np

# ---------------------------------------------------------------------------
# Simulation
# ---------------------------------------------------------------------------
def simulate_data(n_panel_dip, n_target_dip, seq_length=2e6, seed=42):
    """Generate biallelic haplotype matrices from msprime."""
    import msprime

    n_total = n_panel_dip + n_target_dip
    ts = msprime.sim_ancestry(
        samples=n_total, sequence_length=seq_length,
        recombination_rate=1.3e-8, population_size=10000, random_seed=seed)
    ts = msprime.sim_mutations(ts, rate=1.4e-8, random_seed=seed + 1)

    # Extract biallelic sites
    positions_bp = []
    positions_cm = []
    genotypes = []
    for var in ts.variants():
        g = var.genotypes.astype(np.int8)
        if len(np.unique(g)) != 2 or np.any((g != 0) & (g != 1)):
            continue
        positions_bp.append(var.site.position)
        positions_cm.append(var.site.position * 1.3e-8 * 100)
        genotypes.append(g)

    G = np.array(genotypes, dtype=bool)  # (n_sites, n_haps)
    n_panel_haps = 2 * n_panel_dip
    panel_snps = G[:, :n_panel_haps]
    target_snps = G[:, n_panel_haps:]

    cm_pos = np.array(positions_cm)

    return panel_snps, target_snps, cm_pos


# ---------------------------------------------------------------------------
# Benchmark functions
# ---------------------------------------------------------------------------
def bench_set_emission(panel_snps, target_hap, mutation_rate):
    """Benchmark set_emission_probabilities."""
    from threads_arg.fwbw import set_emission_probabilities
    query = target_hap[None, :]
    t0 = time.perf_counter()
    e = set_emission_probabilities(panel_snps, query, mutation_rate)
    return time.perf_counter() - t0, e


def bench_fwbw(panel_subset, target_hap, recomb_rates, mutation_rate):
    """Benchmark full fwbw call."""
    from threads_arg.fwbw import fwbw
    query = target_hap[None, :]
    t0 = time.perf_counter()
    posterior = fwbw(panel_subset, query, recomb_rates, mutation_rate)
    return time.perf_counter() - t0, posterior


def bench_sparsify(posterior, matched_samples, num_samples_panel, num_snps):
    """Benchmark _sparsify_posterior logic (extracted)."""
    from scipy.sparse import csr_array

    t0 = time.perf_counter()
    posterior[posterior <= 1 / num_samples_panel] = 0
    row_sums = posterior.sum(axis=1)
    posterior = posterior / row_sums[:, np.newaxis]
    rows, cols_local = np.nonzero(posterior)
    vals = posterior[rows, cols_local]
    cols_global = matched_samples[cols_local]
    sparse = csr_array(
        (vals, (rows, cols_global)),
        shape=(num_snps, num_samples_panel)
    )
    return time.perf_counter() - t0, sparse


def bench_genotype_sum(posteriors_2d, n_repeats=100):
    """Benchmark vectorized vs list-comprehension genotype sum."""
    # Vectorized (current)
    t0 = time.perf_counter()
    for _ in range(n_repeats):
        g1 = posteriors_2d.sum(axis=1)
    t_vec = (time.perf_counter() - t0) / n_repeats

    # List comprehension (original)
    t0 = time.perf_counter()
    for _ in range(n_repeats):
        g2 = np.array([np.sum(asp) for asp in posteriors_2d])
    t_list = (time.perf_counter() - t0) / n_repeats

    assert np.allclose(g1, g2)
    return t_vec, t_list


def bench_posterior_cache(sparse_posteriors, n_snps, n_repeats=3):
    """Benchmark CachedPosteriorSnps rebuild: vstack vs original loop."""
    from scipy.sparse import vstack as sparse_vstack

    n_targets = len(sparse_posteriors)

    # Current: vstack batch conversion
    t0 = time.perf_counter()
    for _ in range(n_repeats):
        for snp_idx in range(n_snps):
            rows = sparse_vstack([p[[snp_idx]] for p in sparse_posteriors])
            tp = rows.toarray()
            rs = tp.sum(axis=1, keepdims=True)
            tp /= rs
    t_vstack = (time.perf_counter() - t0) / n_repeats / n_snps

    # Original: per-target loop
    t0 = time.perf_counter()
    for _ in range(n_repeats):
        for snp_idx in range(n_snps):
            n_panel = sparse_posteriors[0].shape[1]
            tp = np.empty((n_targets, n_panel), dtype=np.float64)
            for i, p in enumerate(sparse_posteriors):
                row = p[[snp_idx], :].toarray()
                tp[i] = row / np.sum(row)
    t_loop = (time.perf_counter() - t0) / n_repeats / n_snps

    return t_vstack, t_loop


def bench_write_format(n_samples, n_repeats=200):
    """Benchmark VCF GT:DS string formatting: vectorized round vs per-element."""
    rng = np.random.default_rng(42)
    genotypes = rng.random(2 * n_samples)

    # Current: vectorized rint
    t0 = time.perf_counter()
    for _ in range(n_repeats):
        haps1 = genotypes[::2]
        haps2 = genotypes[1::2]
        dosages = haps1 + haps2
        gt1 = np.rint(haps1).astype(int)
        gt2 = np.rint(haps2).astype(int)
        _ = [f"{g1}|{g2}:{d:.3f}".rstrip("0").rstrip(".")
             for g1, g2, d in zip(gt1, gt2, dosages)]
    t_vec = (time.perf_counter() - t0) / n_repeats

    # Original: per-element np.round
    t0 = time.perf_counter()
    for _ in range(n_repeats):
        haps1 = genotypes[::2]
        haps2 = genotypes[1::2]
        dosages = haps1 + haps2
        _ = [f"{np.round(h1):.0f}|{np.round(h2):.0f}:{d:.3f}".rstrip("0").rstrip(".")
             for h1, h2, d in zip(haps1, haps2, dosages)]
    t_orig = (time.perf_counter() - t0) / n_repeats

    return t_vec, t_orig


PRESETS = {
    "small":  {"n_panel": 500,  "n_target": 50,  "seq_length": 2e6, "cond_size": 40},
    "medium": {"n_panel": 1000, "n_target": 100, "seq_length": 5e6, "cond_size": 60},
    "large":  {"n_panel": 2000, "n_target": 200, "seq_length": 10e6, "cond_size": 80},
}


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------
def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--preset", choices=PRESETS.keys(), default=None,
                        help="Use a preset configuration (small/medium/large)")
    parser.add_argument("--n-panel", type=int, default=None,
                        help="Panel diploid count (default 500)")
    parser.add_argument("--n-target", type=int, default=None,
                        help="Target diploid count (default 50)")
    parser.add_argument("--seq-length", type=float, default=None,
                        help="Sequence length in bp (default 2Mb)")
    parser.add_argument("--cond-size", type=int, default=None,
                        help="Conditioning set size per fwbw call (default 40)")
    args = parser.parse_args()

    # Apply preset defaults, then CLI overrides
    preset = PRESETS.get(args.preset, PRESETS["small"])
    n_panel_dip = args.n_panel or preset["n_panel"]
    n_target_dip = args.n_target or preset["n_target"]
    seq_length = args.seq_length or preset["seq_length"]
    cond_size_cfg = args.cond_size or preset["cond_size"]
    n_panel_haps = 2 * n_panel_dip
    n_target_haps = 2 * n_target_dip
    cond_size = min(cond_size_cfg, n_panel_haps)

    print(f"{'=' * 65}")
    print(f"  Imputation Scaling Benchmark")
    print(f"  Panel: {n_panel_dip} dip ({n_panel_haps} hap)")
    print(f"  Target: {n_target_dip} dip ({n_target_haps} hap)")
    print(f"  Seq length: {seq_length/1e6:.1f} Mb")
    print(f"  Conditioning set: {cond_size}")
    print(f"{'=' * 65}")

    # Simulate
    print("\n  Simulating...", end="", flush=True)
    t0 = time.perf_counter()
    panel_snps, target_snps, cm_pos = simulate_data(
        n_panel_dip, n_target_dip, seq_length)
    t_sim = time.perf_counter() - t0
    m, n = panel_snps.shape
    print(f" {t_sim:.1f}s  ({m} sites, {n} panel haps, "
          f"{target_snps.shape[1]} target haps)")

    # Recombination rates
    cm_sizes = np.diff(cm_pos, append=cm_pos[-1] - cm_pos[-2] + cm_pos[-1])
    cm_sizes = np.maximum(cm_sizes, 1e-10)
    Ne = 20_000
    recomb_rates = 1 - np.exp(-4 * Ne * 0.01 * cm_sizes / n_panel_haps)
    mutation_rate = 0.0001

    # Warmup numba
    print("  Warming up numba JIT...", end="", flush=True)
    rng = np.random.default_rng(0)
    cond_idx = rng.choice(n_panel_haps, size=cond_size, replace=False)
    cond_idx.sort()
    h0 = target_snps[:, 0]
    _ = bench_fwbw(panel_snps[:, cond_idx], h0, recomb_rates, mutation_rate)
    print(" done")

    # --- Benchmark set_emission_probabilities ---
    n_emit = min(5, n_target_haps)
    t_emit_total = 0
    for i in range(n_emit):
        t_e, _ = bench_set_emission(panel_snps[:, cond_idx],
                                     target_snps[:, i], mutation_rate)
        t_emit_total += t_e
    t_emit_avg = t_emit_total / n_emit

    # --- Benchmark fwbw ---
    n_fwbw = min(5, n_target_haps)
    t_fwbw_total = 0
    for i in range(n_fwbw):
        ci = rng.choice(n_panel_haps, size=cond_size, replace=False)
        ci.sort()
        t_f, posterior = bench_fwbw(panel_snps[:, ci], target_snps[:, i],
                                     recomb_rates, mutation_rate)
        t_fwbw_total += t_f
    t_fwbw_avg = t_fwbw_total / n_fwbw

    # --- Benchmark sparsify ---
    n_sp = min(5, n_target_haps)
    t_sp_total = 0
    for i in range(n_sp):
        ci = rng.choice(n_panel_haps, size=cond_size, replace=False)
        ci.sort()
        _, post = bench_fwbw(panel_snps[:, ci], target_snps[:, i],
                              recomb_rates, mutation_rate)
        t_s, _ = bench_sparsify(post.copy(), ci, n_panel_haps, m)
        t_sp_total += t_s
    t_sp_avg = t_sp_total / n_sp

    # --- Benchmark genotype sum ---
    dummy_posteriors = rng.random((n_target_haps, cond_size))
    t_vec, t_list = bench_genotype_sum(dummy_posteriors)

    # --- Benchmark posterior cache rebuild ---
    # Build a small set of sparse posteriors for the cache benchmark
    print("  Benchmarking posterior cache...", end="", flush=True)
    n_cache_targets = min(n_target_haps, 20)  # match typical batch size
    sparse_posteriors = []
    for i in range(n_cache_targets):
        ci = rng.choice(n_panel_haps, size=cond_size, replace=False)
        ci.sort()
        _, post = bench_fwbw(panel_snps[:, ci], target_snps[:, i % target_snps.shape[1]],
                              recomb_rates, mutation_rate)
        # Sparsify
        _, sp = bench_sparsify(post.copy(), ci, n_panel_haps, m)
        sparse_posteriors.append(sp)
    n_cache_snps = min(20, m)  # test a subset of SNPs
    t_cache_vstack, t_cache_loop = bench_posterior_cache(sparse_posteriors, n_cache_snps)
    print(" done")

    # --- Benchmark write formatting ---
    n_write_samples = n_target_dip  # diploid target samples
    t_write_vec, t_write_orig = bench_write_format(n_write_samples)

    # --- Projected totals ---
    # In real imputation: n_target_haps fwbw calls, ~m variants with genotype sum
    n_variants = m
    t_fwbw_projected = t_fwbw_avg * n_target_haps
    t_emit_projected = t_emit_avg * n_target_haps
    t_sp_projected = t_sp_avg * n_target_haps
    t_geno_vec_projected = t_vec * n_variants
    t_geno_list_projected = t_list * n_variants
    t_cache_vstack_projected = t_cache_vstack * n_variants
    t_cache_loop_projected = t_cache_loop * n_variants
    t_write_vec_projected = t_write_vec * n_variants
    t_write_orig_projected = t_write_orig * n_variants

    # Print results
    print(f"\n  {'Operation':35s}  {'Per-call':>10s}  {'Projected':>10s}")
    print(f"  {'-'*35}  {'-'*10}  {'-'*10}")
    print(f"  {'set_emission_probabilities':35s}  {t_emit_avg*1000:9.3f}ms  "
          f"{t_emit_projected:9.3f}s")
    print(f"  {'fwbw (forward-backward)':35s}  {t_fwbw_avg*1000:9.3f}ms  "
          f"{t_fwbw_projected:9.3f}s")
    print(f"  {'sparsify_posterior':35s}  {t_sp_avg*1000:9.3f}ms  "
          f"{t_sp_projected:9.3f}s")
    print(f"  {'genotype_sum (vectorized)':35s}  {t_vec*1e6:9.1f}us  "
          f"{t_geno_vec_projected:9.3f}s")
    print(f"  {'genotype_sum (list comp)':35s}  {t_list*1e6:9.1f}us  "
          f"{t_geno_list_projected:9.3f}s")
    print(f"  {'posterior_cache (vstack)':35s}  {t_cache_vstack*1000:9.3f}ms  "
          f"{t_cache_vstack_projected:9.3f}s")
    print(f"  {'posterior_cache (loop)':35s}  {t_cache_loop*1000:9.3f}ms  "
          f"{t_cache_loop_projected:9.3f}s")
    print(f"  {'write_format (vectorized)':35s}  {t_write_vec*1e6:9.1f}us  "
          f"{t_write_vec_projected:9.3f}s")
    print(f"  {'write_format (per-element)':35s}  {t_write_orig*1e6:9.1f}us  "
          f"{t_write_orig_projected:9.3f}s")

    print(f"\n  Projected total for {n_target_haps} targets, {n_variants} variants:")
    t_total_opt = t_fwbw_projected + t_sp_projected + t_geno_vec_projected + t_cache_vstack_projected + t_write_vec_projected
    t_total_orig = t_fwbw_projected + t_sp_projected + t_geno_list_projected + t_cache_loop_projected + t_write_orig_projected
    print(f"    optimized: {t_total_opt:.3f}s")
    print(f"    original:  {t_total_orig:.3f}s")
    print(f"    speedup:   {t_total_orig/t_total_opt:.2f}x")
    print()


if __name__ == "__main__":
    main()
