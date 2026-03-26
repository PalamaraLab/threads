"""Benchmark tree multiply on cached .threads files from threads-inferred msprime ARGs.

Loads pre-computed .threads files from bench_data/ (created by simulate_threads.py).
Reports prepare, left multiply, right multiply timings and scaling exponents.
"""
import os
import sys
import time
import glob
import numpy as np
from threads_arg.serialization import load_instructions

SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
CACHE_DIR = os.path.join(SCRIPT_DIR, "bench_data")


def bench_one(path):
    """Benchmark a single .threads file."""
    ti = load_instructions(path)
    n = ti.num_samples
    m = ti.num_sites

    # Prepare
    t0 = time.perf_counter()
    ti.prepare_tree_multiply()
    t_prep = time.perf_counter() - t0

    # Vectors
    rng = np.random.default_rng(42)
    xl = rng.standard_normal(n)
    xr = rng.standard_normal(m)

    # Warm up
    ti.left_multiply_tree(xl)
    ti.right_multiply_tree(xr)

    # Calibrate reps
    t0 = time.perf_counter()
    ti.left_multiply_tree(xl)
    t_single = time.perf_counter() - t0
    reps = max(1, min(200, int(2.0 / max(t_single, 0.0001))))

    t0 = time.perf_counter()
    for _ in range(reps):
        ti.left_multiply_tree(xl)
    t_left = (time.perf_counter() - t0) / reps

    t0 = time.perf_counter()
    for _ in range(reps):
        ti.right_multiply_tree(xr)
    t_right = (time.perf_counter() - t0) / reps

    return {"n": n, "m": m, "prep_s": t_prep, "left_ms": t_left * 1000, "right_ms": t_right * 1000}


def correctness_check(path):
    """Verify tree multiply matches dense multiply."""
    ti = load_instructions(path)
    n = ti.num_samples
    m = ti.num_sites
    rng = np.random.default_rng(99)
    xl = rng.standard_normal(n)
    xr = rng.standard_normal(m)

    dl = np.array(ti.left_multiply(xl.tolist()))
    tl = np.array(ti.left_multiply_tree(xl))
    dr = np.array(ti.right_multiply(xr.tolist()))
    tr = np.array(ti.right_multiply_tree(xr))

    el = np.max(np.abs(dl - tl)) if len(dl) > 0 else 0
    er = np.max(np.abs(dr - tr)) if len(dr) > 0 else 0
    ok = el < 1e-8 and er < 1e-8
    print(f"  n={n:5d} m={m:5d} left_err={el:.2e} right_err={er:.2e} {'PASS' if ok else 'FAIL'}")
    return ok


def main():
    files = sorted(glob.glob(os.path.join(CACHE_DIR, "n*.threads")),
                   key=lambda f: int(os.path.basename(f).replace("n", "").replace(".threads", "")))
    if not files:
        print(f"No .threads files in {CACHE_DIR}. Run simulate_threads.py first.")
        return

    print(f"Found {len(files)} cached .threads files\n")

    # Correctness (only for small n where dense multiply is feasible)
    print("=== Correctness ===")
    for f in files:
        ti_tmp = load_instructions(f)
        if ti_tmp.num_samples > 5000:
            print(f"  n={ti_tmp.num_samples:5d} skipped (dense too slow)")
            continue
        if not correctness_check(f):
            print("CORRECTNESS FAILURE")
            return

    # Scaling
    print(f"\n=== Scaling (threads-inferred, 1Mb region) ===")
    header = f"{'n':>8} {'m':>8} {'prep_s':>10} {'left_ms':>10} {'right_ms':>10}"
    print(header)
    print("-" * len(header))

    results = []
    for f in files:
        r = bench_one(f)
        results.append(r)
        print(f"{r['n']:8d} {r['m']:8d} {r['prep_s']:10.3f} {r['left_ms']:10.2f} {r['right_ms']:10.2f}")

    # Scaling exponents from last 2 points
    if len(results) >= 2:
        print("\n=== Scaling exponents ===")
        r2, r3 = results[-2], results[-1]
        log_r = np.log(r3["n"] / r2["n"])
        for key, label in [("prep_s", "prepare"), ("left_ms", "left"), ("right_ms", "right")]:
            if r2[key] > 0 and r3[key] > 0:
                alpha = np.log(r3[key] / r2[key]) / log_r
                ext_val = r3[key] * (500000 / r3["n"]) ** alpha
                if key == "prep_s":
                    print(f"  {label:>8s} ~ n^{alpha:.2f}  extrapolated 500K: {ext_val:.1f}s")
                else:
                    print(f"  {label:>8s} ~ n^{alpha:.2f}  extrapolated 500K: {ext_val/1000:.1f}s")


if __name__ == "__main__":
    main()
