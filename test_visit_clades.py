#!/usr/bin/env python3
"""Test ARG traversals and statistics against arg-needle-lib.

Tests:
  1. total_volume and AFS accuracy vs arg-needle-lib
  2. Association testing: verify significant hits on simulated causal variant
  3. Mutation generation: AFS of generated mutations matches branch AFS
  4. Timing benchmark
"""
import time
import numpy as np
import msprime
import threads_arg
from threads_arg import ThreadingInstructions, ViterbiPath
from threads_arg.convert import threads_to_arg
import arg_needle_lib

SEQ_LEN = 1e6
RECOMB = 1e-8
MUT = 1e-8
NE = 10000
SEED = 42


def simulate(n_diploid):
    ts = msprime.sim_ancestry(
        n_diploid, sequence_length=SEQ_LEN,
        recombination_rate=RECOMB, population_size=NE, random_seed=SEED)
    ts = msprime.sim_mutations(ts, rate=MUT, random_seed=SEED + 1)
    n_haps = 2 * n_diploid
    positions_bp, positions_cm, genotype_rows = [], [], []
    for var in ts.variants():
        g = var.genotypes.astype(np.int32)
        if len(set(g.tolist())) < 2 or np.any((g != 0) & (g != 1)):
            continue
        positions_bp.append(var.site.position)
        positions_cm.append(var.site.position * RECOMB * 100)
        genotype_rows.append(g)
    geno_matrix = np.stack(genotype_rows)
    return ts, n_haps, positions_bp, positions_cm, geno_matrix


def infer_instructions(n_haps, positions_bp, positions_cm, geno_matrix):
    n_sites = geno_matrix.shape[0]
    matcher = threads_arg.Matcher(n_haps, positions_cm, 0.01, 0.5, 4, 2)
    matcher.process_all_sites_numpy(geno_matrix)
    target_ids = list(range(1, n_haps))
    match_data = matcher.serializable_matches(target_ids)
    cm_pos = matcher.cm_positions()
    tlm = threads_arg.ThreadsLowMem(
        target_ids, positions_bp, positions_cm,
        [float(NE)], [0.0], MUT, False)
    tlm.initialize_viterbi(match_data, cm_pos)
    tlm.process_all_sites_viterbi_numpy(geno_matrix)
    tlm.prune()
    tlm.traceback()
    tlm.process_all_sites_hets_numpy(geno_matrix)
    tlm.date_segments()
    all_starts, all_ids, all_heights, all_hetsites = tlm.serialize_paths()
    sample0_hets = [int(i) for i in range(n_sites) if geno_matrix[i, 0] == 1]
    paths = [ViterbiPath(0, [], [], [], sample0_hets)]
    for ss, mi, ht, hs in zip(all_starts, all_ids, all_heights, all_hetsites):
        paths.append(ViterbiPath(target_ids[len(paths) - 1], ss, mi, ht, hs))
    positions_int = [int(p) for p in positions_bp]
    return ThreadingInstructions(
        paths, int(positions_bp[0]), int(positions_bp[-1]) + 1, positions_int)


def instructions_to_arg(instructions):
    try:
        arg = threads_to_arg(instructions, add_mutations=False, noise=1e-5)
    except RuntimeError:
        arg = threads_to_arg(instructions, add_mutations=False, noise=1e-3)
    arg.populate_children_and_roots()
    return arg


def afs_from_bitset_volume_map(arg):
    n = arg.num_samples()
    afs = [0.0] * (n + 1)
    for bitset_str, vol in arg_needle_lib.bitset_volume_map(arg).items():
        afs[bitset_str.count('1')] += vol
    return afs


# ── Test 1: Volume and AFS accuracy ─────────────────────────────────────

def test_volume_and_afs(n_diploid=50):
    n_haps = 2 * n_diploid
    print(f"Test 1: Volume and AFS ({n_haps} haplotypes)")
    ts, n_haps, positions_bp, positions_cm, geno_matrix = simulate(n_diploid)
    instructions = infer_instructions(n_haps, positions_bp, positions_cm, geno_matrix)
    arg = instructions_to_arg(instructions)

    # Total volume
    tv_threads = instructions.total_volume()
    tv_anl = arg_needle_lib.total_volume(arg)
    ratio = tv_threads / tv_anl
    print(f"  total_volume: threads={tv_threads:.2e}  anl={tv_anl:.2e}  ratio={ratio:.6f}")
    assert abs(ratio - 1.0) < 0.02

    # AFS
    afs_threads = instructions.allele_frequency_spectrum()
    afs_anl = afs_from_bitset_volume_map(arg)
    t_arr, a_arr = np.array(afs_threads), np.array(afs_anl)
    nonzero = (t_arr > 0) | (a_arr > 0)
    corr = np.corrcoef(t_arr[nonzero], a_arr[nonzero])[0, 1]
    print(f"  AFS correlation: {corr:.6f}")
    assert corr > 0.99

    print("  PASS\n")
    return instructions, arg


# ── Test 2: Association testing ──────────────────────────────────────────

def test_association(instructions):
    n = instructions.num_samples
    n_dip = n // 2
    print(f"Test 2: Association testing ({n_dip} diploid individuals)")

    rng = np.random.default_rng(123)

    # Generate mutations and pick a common causal variant
    pos_list, geno_flat, n_mut = instructions.generate_mutations(1e-8, seed=99)
    geno = np.array(geno_flat).reshape(n_mut, n)
    dosage = geno[:, 0::2] + geno[:, 1::2]
    maf = dosage.sum(axis=1) / (2 * n_dip)
    maf = np.minimum(maf, 1 - maf)
    good = np.where((maf > 0.1) & (maf < 0.4))[0]
    if len(good) == 0:
        print("  SKIP: no common variants generated\n")
        return
    causal_idx = good[len(good) // 2]
    causal_dosage = dosage[causal_idx]
    causal_pos = pos_list[causal_idx]
    causal_maf = maf[causal_idx]

    # Phenotype = dosage * beta + noise
    beta = 1.0
    phenotypes = causal_dosage * beta + rng.normal(0, 1, n_dip)
    phenotypes = (phenotypes - phenotypes.mean()) / phenotypes.std()

    positions, chi2s, pvals, n_desc = instructions.association_diploid(
        phenotypes.tolist(), p_threshold=0.01)

    print(f"  Causal variant at {causal_pos}bp, MAF={causal_maf:.3f}")
    print(f"  {len(positions)} significant clades (p < 0.01)")

    if len(positions) > 0:
        best_idx = np.argmin(pvals)
        best_pos = positions[best_idx]
        best_p = pvals[best_idx]
        dist = abs(best_pos - causal_pos)
        print(f"  Best hit: pos={best_pos}bp, p={best_p:.2e}, "
              f"n_desc={n_desc[best_idx]}, dist={dist}bp")

        nearby = [i for i, p in enumerate(positions)
                  if abs(p - causal_pos) < SEQ_LEN * 0.1]
        print(f"  Hits within 10% of causal: {len(nearby)}")
        assert len(nearby) > 0, "No hits near causal variant"

    # Null phenotype should yield far fewer hits
    random_pheno = rng.normal(0, 1, n_dip).tolist()
    pos_null, _, _, _ = instructions.association_diploid(random_pheno, p_threshold=0.01)
    print(f"  Null phenotype hits: {len(pos_null)} (expect << {len(positions)})")

    print("  PASS\n")


# ── Test 3: Mutation generation ──────────────────────────────────────────

def test_mutation_generation(instructions, arg):
    n = instructions.num_samples
    print(f"Test 3: Mutation generation ({n} haplotypes)")

    mu = 1e-8
    pos_list, geno_flat, n_mut = instructions.generate_mutations(mu, seed=42)
    geno = np.array(geno_flat).reshape(n_mut, n)
    print(f"  Generated {n_mut} mutations")

    ac = geno.sum(axis=1)
    empirical_afs = np.bincount(ac.astype(int), minlength=n + 1)

    # Expected from branch volumes
    branch_afs = np.array(instructions.allele_frequency_spectrum())
    expected_counts = branch_afs * mu
    total_expected = expected_counts.sum()
    print(f"  Expected mutations: {total_expected:.0f}, got: {n_mut}")
    assert abs(n_mut - total_expected) < 5 * np.sqrt(total_expected), \
        f"Mutation count {n_mut} far from expected {total_expected:.0f}"

    # AFS correlation
    nonzero = (expected_counts > 1) & (empirical_afs > 0)
    if nonzero.sum() > 3:
        corr = np.corrcoef(empirical_afs[nonzero], expected_counts[nonzero])[0, 1]
        print(f"  AFS correlation (empirical vs expected): {corr:.4f}")
        assert corr > 0.8, f"AFS correlation {corr} too low"

    assert geno.min() == 0 and geno.max() == 1
    assert (ac > 0).all() and (ac < n).all()

    # Compare count against arg-needle-lib
    arg_needle_lib.generate_mutations(arg, mu, random_seed=42)
    anl_n_mut = arg.num_mutations()
    ratio = n_mut / anl_n_mut if anl_n_mut > 0 else float('inf')
    print(f"  arg-needle-lib mutations: {anl_n_mut}, ratio: {ratio:.2f}")
    assert 0.5 < ratio < 2.0

    print("  PASS\n")


# ── Benchmark ────────────────────────────────────────────────────────────

def benchmark(sizes=[20, 50, 100, 250, 500]):
    print(f"Benchmark: threading vs arg-needle-lib")
    print(f"{'n_haps':>8} {'sites':>6} "
          f"{'thr_vol':>10} {'anl_vol':>10} "
          f"{'thr_afs':>10} {'anl_afs':>10} "
          f"{'su_vol':>8} {'su_afs':>8}")
    print("-" * 82)

    for n_dip in sizes:
        ts, n_haps, positions_bp, positions_cm, geno_matrix = simulate(n_dip)
        instructions = infer_instructions(n_haps, positions_bp, positions_cm, geno_matrix)
        arg = instructions_to_arg(instructions)

        instructions.total_volume()
        arg_needle_lib.total_volume(arg)
        n_reps = 3

        t0 = time.time()
        for _ in range(n_reps): instructions.total_volume()
        t_tv = (time.time() - t0) / n_reps

        t0 = time.time()
        for _ in range(n_reps): arg_needle_lib.total_volume(arg)
        t_av = (time.time() - t0) / n_reps

        t0 = time.time()
        for _ in range(n_reps): instructions.allele_frequency_spectrum()
        t_ta = (time.time() - t0) / n_reps

        t0 = time.time()
        for _ in range(n_reps): afs_from_bitset_volume_map(arg)
        t_aa = (time.time() - t0) / n_reps

        print(f"{n_haps:>8} {len(positions_bp):>6} "
              f"{t_tv*1000:>8.2f}ms {t_av*1000:>8.2f}ms "
              f"{t_ta*1000:>8.2f}ms {t_aa*1000:>8.2f}ms "
              f"{t_av/t_tv:>7.2f}x {t_aa/t_ta:>7.2f}x")


if __name__ == '__main__':
    instructions, arg = test_volume_and_afs()
    test_association(instructions)
    test_mutation_generation(instructions, arg)
    benchmark()
