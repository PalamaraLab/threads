#!/usr/bin/env python3
"""
Quick compression & multiply benchmark: threads vs GRG vs RePair.
Uses pre-cached objects only — no rebuilding.

Usage:
  python test/bench_compression_multiply.py [--sizes 100,1000,10000] [--reps 5]
  python test/bench_compression_multiply.py --dc --sizes 100  # also run data consistency
"""

import os, sys, time, argparse, json, tempfile
import numpy as np

sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..', 'src'))
_base = os.path.join(os.path.dirname(__file__), '..', '..')
sys.path.insert(0, os.path.join(_base, 'genrepair'))
sys.path.insert(0, os.path.join(_base, 'other', 'grgl'))

from threads_arg.serialization import load_instructions, serialize_instructions

TC = os.path.join(os.path.dirname(os.path.abspath(__file__)), "threads_cache")
DC = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "threads_cache")


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
    ap.add_argument("--sizes", default="100,1000,10000")
    ap.add_argument("--reps", type=int, default=5)
    ap.add_argument("--dc", action="store_true", help="Run data consistency (slow)")
    args = ap.parse_args()
    sizes = [int(s) for s in args.sizes.split(",")]

    try:
        import pygrgl; have_grg = True
    except ImportError:
        have_grg = False; print("no pygrgl")

    try:
        import genrepair; have_rpr = True
    except ImportError:
        have_rpr = False; print("no genrepair")

    results = []
    for n in sizes:
        # Find files
        tp = os.path.join(TC, f"sim_{n}dip.threads")
        if not os.path.exists(tp):
            tp = os.path.join(DC, f"sim_{n}.threads")
        if not os.path.exists(tp):
            print(f"{n}: skip"); continue

        inst = load_instructions(tp)
        ns, nm = inst.num_samples, inst.num_sites
        print(f"\n--- {n} dip ({ns} hap, {nm} sites) ---")

        # Sizes
        dense = nm * 2 * n
        thr_sz = os.path.getsize(tp)
        grg_p = os.path.join(TC, f"grg_{n}dip.grg")
        grg_sz = os.path.getsize(grg_p) if os.path.exists(grg_p) else None
        rpr_p = os.path.join(DC, f"sim_{n}.repair")
        rpr_sz = None
        if have_rpr:
            for ext in [".val", ".vc.C.ansf.1", ".vc.R.iv"]:
                if not os.path.exists(rpr_p + ext):
                    rpr_sz = None; break
                rpr_sz = (rpr_sz or 0) + os.path.getsize(rpr_p + ext)

        print(f"  size: dense={sz(dense)}  threads={sz(thr_sz)}({dense/thr_sz:.0f}x)  "
              f"grg={sz(grg_sz)}({dense/grg_sz:.0f}x)" if grg_sz else f"  size: dense={sz(dense)}  threads={sz(thr_sz)}({dense/thr_sz:.0f}x)",
              end="")
        if rpr_sz:
            print(f"  repair={sz(rpr_sz)}({dense/rpr_sz:.0f}x)", end="")
        print()

        # Multiply
        rng = np.random.default_rng(42)
        xs = rng.normal(0, 1, nm).tolist()
        xh = rng.normal(0, 1, ns).tolist()

        inst.prepare_tree_multiply()
        tr = bench(lambda: inst.right_multiply_tree(xs), args.reps)
        tl = bench(lambda: inst.left_multiply_tree(xh), args.reps)

        inst.materialize_genotypes()
        dr = bench(lambda: inst.right_multiply(xs), args.reps)
        dl = bench(lambda: inst.left_multiply(xh), args.reps)

        r = dict(n=n, ns=ns, nm=nm, dense_bytes=dense, thr_bytes=thr_sz,
                 grg_bytes=grg_sz, rpr_bytes=rpr_sz,
                 tree_R=tr, tree_L=tl, dense_R=dr, dense_L=dl)

        print(f"  right: tree={tr:.2f}ms  dense={dr:.2f}ms({dr/tr:.1f}x)", end="")

        if have_grg and os.path.exists(grg_p):
            grg = pygrgl.load_immutable_grg(grg_p)
            xg = rng.normal(0, 1, grg.num_mutations).astype(np.float32)
            gr = bench(lambda: pygrgl.dot_product(grg, xg, pygrgl.DOWN), args.reps)
            xgu = rng.normal(0, 1, grg.num_samples).astype(np.float32)
            gl = bench(lambda: pygrgl.dot_product(grg, xgu, pygrgl.UP), args.reps)
            r['grg_R'] = gr; r['grg_L'] = gl
            print(f"  grg={gr:.2f}ms({gr/tr:.1f}x)", end="")

        if have_rpr and rpr_sz:
            cm = genrepair.CompressedMatrix(rpr_p, nm, 2*n, False, False)
            xr = rng.normal(0, 1, cm.rows).astype(np.float32)
            xl = rng.normal(0, 1, cm.cols).astype(np.float32)
            rr = bench(lambda: cm.left_multiply(xr), args.reps)
            rl = bench(lambda: cm.multiply(xl), args.reps)
            r['rpr_R'] = rr; r['rpr_L'] = rl
            print(f"  repair={rr:.2f}ms({rr/tr:.1f}x)", end="")
        print()

        print(f"  left:  tree={tl:.2f}ms  dense={dl:.2f}ms({dl/tl:.1f}x)", end="")
        if 'grg_L' in r:
            print(f"  grg={r['grg_L']:.2f}ms({r['grg_L']/tl:.1f}x)", end="")
        if 'rpr_L' in r:
            print(f"  repair={r['rpr_L']:.2f}ms({r['rpr_L']/tl:.1f}x)", end="")
        print()

        # Data consistency (optional, slow)
        if args.dc:
            gp = os.path.join(DC, f"sim_{n}.genotypes.npz")
            if os.path.exists(gp):
                from threads_arg import AgeEstimator, GenotypeIterator, ConsistencyWrapper
                print(f"  DC: building...", end="", flush=True)
                inst_dc = load_instructions(tp)
                t0 = time.perf_counter()
                ae = AgeEstimator(inst_dc)
                gi = GenotypeIterator(inst_dc)
                while gi.has_next_genotype():
                    ae.process_site(np.array(gi.next_genotype()))
                ages = ae.get_inferred_ages()
                cw = ConsistencyWrapper(inst_dc, ages)
                gi2 = GenotypeIterator(inst_dc)
                while gi2.has_next_genotype():
                    cw.process_site(gi2.next_genotype())
                dc_inst = cw.get_consistent_instructions()
                dc_t = time.perf_counter() - t0
                with tempfile.NamedTemporaryFile(suffix='.threads', delete=False) as f:
                    dcp = f.name
                serialize_instructions(dc_inst, dcp, allele_ages=ages)
                dc_sz = os.path.getsize(dcp)
                dc_inst2 = load_instructions(dcp)
                dc_inst2.prepare_tree_multiply()
                xsd = rng.normal(0, 1, nm).tolist()
                xhd = rng.normal(0, 1, dc_inst2.num_samples).tolist()
                dcr = bench(lambda: dc_inst2.right_multiply_tree(xsd), args.reps)
                dcl = bench(lambda: dc_inst2.left_multiply_tree(xhd), args.reps)
                os.unlink(dcp)
                r['dc_bytes'] = dc_sz; r['dc_build_s'] = dc_t
                r['dc_R'] = dcr; r['dc_L'] = dcl
                print(f" {dc_t:.1f}s  size={sz(dc_sz)}({dc_sz/thr_sz:.2f}x std)  "
                      f"right={dcr:.2f}ms({dcr/tr:.1f}x)  left={dcl:.2f}ms({dcl/tl:.1f}x)")

        results.append(r)

    out = os.path.join(DC, "bench_compression_multiply.json")
    with open(out, 'w') as f:
        json.dump(results, f, indent=2)
    print(f"\nSaved to {out}")


if __name__ == "__main__":
    main()
