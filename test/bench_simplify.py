#!/usr/bin/env python3
"""
Benchmark: RLE vs tree vs dense multiply from .threads format.

Usage:
  python test/bench_simplify.py [--sizes 100,1000,10000] [--reps 5]
"""

import os, sys, time, argparse, json, tempfile
import numpy as np

sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..', 'src'))

from threads_arg.serialization import load_instructions, serialize_instructions

TC = os.path.join(os.path.dirname(os.path.abspath(__file__)), "threads_cache")


def bench(fn, reps=5):
    fn()
    times = []
    for _ in range(reps):
        t0 = time.perf_counter()
        fn()
        times.append((time.perf_counter() - t0) * 1000)
    return np.median(times)


def sz(b):
    if b is None: return "—"
    if b < 1024: return f"{b}B"
    if b < 1024**2: return f"{b/1024:.0f}K"
    return f"{b/1024/1024:.1f}M"


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--sizes", default="50,100,250,500,1000,2000,5000,10000")
    ap.add_argument("--reps", type=int, default=5)
    args = ap.parse_args()
    sizes = [int(s) for s in args.sizes.split(",")]

    try:
        import pygrgl; have_grg = True
    except ImportError:
        have_grg = False; print("(no pygrgl)")

    results = []
    for n in sizes:
        tp = os.path.join(TC, f"sim_{n}dip.threads")
        if not os.path.exists(tp):
            print(f"{n}: skip"); continue

        inst = load_instructions(tp)
        ns, nm = inst.num_samples, inst.num_sites
        orig_segs = sum(len(s) for s in inst.all_starts())
        orig_mm = sum(len(m) for m in inst.all_mismatches())
        orig_sz = os.path.getsize(tp)

        # Simplification file size
        simp = inst.simplify_for_multiply()
        simp_segs = sum(len(s) for s in simp.all_starts())
        with tempfile.NamedTemporaryFile(suffix='.threads', delete=False) as f:
            simp_path = f.name
        serialize_instructions(simp, simp_path)
        simp_sz = os.path.getsize(simp_path)
        os.unlink(simp_path)

        dense_sz = nm * ns  # 1 byte per entry
        sz_reduction = 100 * (1 - simp_sz / orig_sz)

        print(f"\n--- {n} dip ({ns} hap, {nm} sites, {orig_segs} segs, {orig_mm} mm) ---")
        print(f"  size: threads={sz(orig_sz)}  simplified={sz(simp_sz)} ({sz_reduction:.0f}% smaller)  dense={sz(dense_sz)}")

        rng = np.random.default_rng(42)
        xs = rng.normal(0, 1, nm).tolist()
        xh = rng.normal(0, 1, ns).tolist()

        # Prepare times
        inst_rle = load_instructions(tp)
        t0 = time.perf_counter()
        inst_rle.prepare_rle_multiply()
        rle_prep = (time.perf_counter() - t0) * 1000

        inst_tree = load_instructions(tp)
        t0 = time.perf_counter()
        inst_tree.prepare_tree_multiply()
        tree_prep = (time.perf_counter() - t0) * 1000

        print(f"  prepare: rle={rle_prep:.1f}ms  tree={tree_prep:.1f}ms")

        # Multiply benchmarks
        rle_r = bench(lambda: inst_rle.right_multiply_rle(xs), args.reps)
        rle_l = bench(lambda: inst_rle.left_multiply_rle(xh), args.reps)
        tree_r = bench(lambda: inst_tree.right_multiply_tree(xs), args.reps)
        tree_l = bench(lambda: inst_tree.left_multiply_tree(xh), args.reps)

        inst.materialize_genotypes()
        dense_r = bench(lambda: inst.right_multiply(xs), args.reps)
        dense_l = bench(lambda: inst.left_multiply(xh), args.reps)

        r = dict(
            n=n, ns=ns, nm=nm,
            segs=orig_segs, simp_segs=simp_segs, mm=orig_mm,
            orig_bytes=orig_sz, simp_bytes=simp_sz, dense_bytes=dense_sz,
            rle_prep=rle_prep, tree_prep=tree_prep,
            rle_R=rle_r, rle_L=rle_l,
            tree_R=tree_r, tree_L=tree_l,
            dense_R=dense_r, dense_L=dense_l,
        )

        print(f"  right:  rle={rle_r:.2f}  tree={tree_r:.2f}  dense={dense_r:.2f}ms")
        print(f"  left:   rle={rle_l:.2f}  tree={tree_l:.2f}  dense={dense_l:.2f}ms")

        # GRG
        grg_p = os.path.join(TC, f"grg_{n}dip.grg")
        if have_grg and os.path.exists(grg_p):
            grg = pygrgl.load_immutable_grg(grg_p)
            grg_sz = os.path.getsize(grg_p)
            xg = rng.normal(0, 1, grg.num_mutations).astype(np.float32)
            xgu = rng.normal(0, 1, grg.num_samples).astype(np.float32)
            gr = bench(lambda: pygrgl.dot_product(grg, xg, pygrgl.DOWN), args.reps)
            gl = bench(lambda: pygrgl.dot_product(grg, xgu, pygrgl.UP), args.reps)
            r['grg_bytes'] = grg_sz; r['grg_R'] = gr; r['grg_L'] = gl
            print(f"  grg:    right={gr:.2f}  left={gl:.2f}ms  size={sz(grg_sz)}")

        results.append(r)

    out = os.path.join(TC, "bench_simplify.json")
    with open(out, 'w') as f:
        json.dump(results, f, indent=2)
    print(f"\nSaved to {out}")

    # Summary table
    print(f"\n{'n':>6} | {'rle_R':>7} {'tree_R':>7} {'dense_R':>7} | {'rle_L':>7} {'tree_L':>7} {'dense_L':>7} | {'rle/tree':>8}")
    print("-" * 80)
    for r in results:
        rt = r['tree_R']/r['rle_R']; lt = r['tree_L']/r['rle_L']
        print(f"{r['n']:>6} | {r['rle_R']:>6.1f}  {r['tree_R']:>6.1f}  {r['dense_R']:>6.1f}  | "
              f"{r['rle_L']:>6.1f}  {r['tree_L']:>6.1f}  {r['dense_L']:>6.1f}  | "
              f"R:{rt:.1f}x L:{lt:.1f}x")


if __name__ == "__main__":
    main()
