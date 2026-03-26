"""Compare committed vs working tree shuttle multiply.

Runs benchmarks with different builds via importlib to bypass site-packages.
"""
import os, sys, time, glob, json
import importlib.util
import numpy as np

SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
CACHE_DIR = os.path.join(SCRIPT_DIR, "bench_data")


def load_build(so_path):
    """Force-load a specific .so file, bypassing site-packages."""
    spec = importlib.util.spec_from_file_location(
        "threads_arg_python_bindings", so_path)
    m = importlib.util.module_from_spec(spec)
    sys.modules["threads_arg_python_bindings"] = m
    spec.loader.exec_module(m)
    return m


def bench_all(files, so_path):
    """Benchmark a specific build on all cached .threads files."""
    import subprocess
    python = sys.executable
    results = []
    for f in files:
        # Run each size in a fresh subprocess for clean RSS measurement
        code = f"""
import importlib.util, sys, time, resource, json
import numpy as np
spec = importlib.util.spec_from_file_location(
    "threads_arg_python_bindings", "{so_path}")
m = importlib.util.module_from_spec(spec)
sys.modules["threads_arg_python_bindings"] = m
spec.loader.exec_module(m)
from threads_arg.serialization import load_instructions
ti = load_instructions("{f}")
n, m_sites = ti.num_samples, ti.num_sites
rng = np.random.default_rng(42)
xl = rng.standard_normal(n)
xr = rng.standard_normal(m_sites)

# Prepare + memory
rss0 = resource.getrusage(resource.RUSAGE_SELF).ru_maxrss
t0 = time.perf_counter()
ti.prepare_tree_multiply()
t_prep = (time.perf_counter() - t0) * 1000
rss1 = resource.getrusage(resource.RUSAGE_SELF).ru_maxrss
mem_prep = max(0, rss1 - rss0) / (1024 * 1024)

# Memory for left multiply (first call may allocate W matrix)
rss2 = resource.getrusage(resource.RUSAGE_SELF).ru_maxrss
ti.left_multiply_tree(xl)
rss3 = resource.getrusage(resource.RUSAGE_SELF).ru_maxrss
mem_left = max(0, rss3 - rss2) / (1024 * 1024)

# Warmup
ti.right_multiply_tree(xr)

# Calibrate
t0 = time.perf_counter()
ti.left_multiply_tree(xl)
t_single = time.perf_counter() - t0
reps = max(1, min(200, int(2.0 / max(t_single, 0.0001))))

t0 = time.perf_counter()
for _ in range(reps): ti.left_multiply_tree(xl)
t_left = (time.perf_counter() - t0) / reps * 1000

t0 = time.perf_counter()
for _ in range(reps): ti.right_multiply_tree(xr)
t_right = (time.perf_counter() - t0) / reps * 1000

print(json.dumps({{"n": n, "m": m_sites, "prep_ms": t_prep, "left_ms": t_left,
                   "right_ms": t_right, "mem_prep": mem_prep, "mem_left": mem_left}}))
"""
        p = subprocess.run([python, "-c", code], capture_output=True, text=True,
                           env={**os.environ, "PYTHONPATH": os.path.join(SCRIPT_DIR, "src")})
        if p.returncode != 0:
            print(f"  Failed for {f}:\n{p.stderr[:500]}")
            return None
        results.append(json.loads(p.stdout.strip().split("\n")[-1]))
    return results


def print_results(label, results):
    print(f"\n=== {label} ===")
    print(f"{'n':>8} {'m':>8} {'prep_ms':>10} {'left_ms':>10} {'right_ms':>10} {'mem_prep':>10} {'mem_left':>10}")
    print("-" * 80)
    for r in results:
        print(f"{r['n']:8d} {r['m']:8d} {r['prep_ms']:10.1f} {r['left_ms']:10.2f} "
              f"{r['right_ms']:10.2f} {r['mem_prep']:10.1f} {r['mem_left']:10.1f}")


def main():
    files = sorted(glob.glob(os.path.join(CACHE_DIR, "n*.threads")),
                   key=lambda f: int(os.path.basename(f).replace("n", "").replace(".threads", "")))
    if not files:
        print(f"No .threads files in {CACHE_DIR}")
        return

    committed_so = os.path.join(SCRIPT_DIR, "build_committed", "src",
                                "threads_arg_python_bindings.cpython-310-darwin.so")
    working_so = os.path.join(SCRIPT_DIR, "build_bench", "src",
                              "threads_arg_python_bindings.cpython-310-darwin.so")

    for so, label in [(committed_so, "committed"), (working_so, "working")]:
        if not os.path.exists(so):
            print(f"Missing {label} build: {so}")
            return

    print("Benchmarking COMMITTED version...")
    r_committed = bench_all(files, committed_so)
    if not r_committed:
        return

    print("Benchmarking WORKING version...")
    r_working = bench_all(files, working_so)
    if not r_working:
        return

    print_results("COMMITTED (visit-clades HEAD)", r_committed)
    print_results("WORKING (uncommitted changes)", r_working)

    print(f"\n=== SPEEDUP (committed / working) ===")
    print(f"{'n':>8} {'prep':>10} {'left':>10} {'right':>10} {'mem_left':>12}")
    print("-" * 56)
    for c, w in zip(r_committed, r_working):
        sp = c["prep_ms"] / w["prep_ms"] if w["prep_ms"] > 0 else 0
        sl = c["left_ms"] / w["left_ms"] if w["left_ms"] > 0 else 0
        sr = c["right_ms"] / w["right_ms"] if w["right_ms"] > 0 else 0
        if w["mem_left"] > 0:
            sm = f"{c['mem_left'] / w['mem_left']:.0f}x"
        elif c["mem_left"] > 0:
            sm = "inf"
        else:
            sm = "~"
        print(f"{c['n']:8d} {sp:9.2f}x {sl:9.2f}x {sr:9.2f}x {sm:>12s}")


if __name__ == "__main__":
    main()
