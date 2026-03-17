#!/usr/bin/env python3
"""
Imputation benchmark: correctness (bit-identical output), speed, and memory.

Usage:
    python test/bench_impute.py                # run baseline, save reference
    python test/bench_impute.py --compare      # run optimized, compare to reference
    python test/bench_impute.py --repeat 3     # best-of-3 runs

The first run (without --compare) records:
  - The full VCF output (reference for bit-identity checks)
  - Wall time for each pipeline stage
  - Peak RSS memory

The second run (with --compare) re-runs imputation and checks:
  1. Output is bit-identical to the reference (ignoring the date header line)
  2. Wall time per stage (prints speedup ratios)
  3. Peak memory (prints reduction)

Results are saved to test/bench_impute_results/ so you can compare across
code changes without re-running the baseline.
"""
import argparse
import json
import logging
import os
import re
import resource
import sys
import tempfile
import time

from pathlib import Path

# ---------------------------------------------------------------------------
# Paths
# ---------------------------------------------------------------------------
BASE_DIR = Path(__file__).parent.parent
TEST_DATA_DIR = BASE_DIR / "test" / "data"
RESULTS_DIR = BASE_DIR / "test" / "bench_impute_results"

# Inputs (pre-generated snapshot fixtures)
PANEL_VCF = TEST_DATA_DIR / "panel.vcf.gz"
TARGET_VCF = TEST_DATA_DIR / "target.vcf.gz"
GMAP = TEST_DATA_DIR / "gmap_04.map"
MUT = TEST_DATA_DIR / "expected_mapping_snapshot.mut"
DEMO = TEST_DATA_DIR / "CEU_unscaled.demo"
REGION = "1:400000-600000"

# Reference snapshot for bit-identity (the existing regression fixture)
EXPECTED_VCF = TEST_DATA_DIR / "expected_impute_snapshot.vcf"


# ---------------------------------------------------------------------------
# Measurement helpers
# ---------------------------------------------------------------------------
def peak_rss_mb():
    """Current peak RSS in MB (macOS returns bytes, Linux returns KB)."""
    ru = resource.getrusage(resource.RUSAGE_SELF)
    if sys.platform == "darwin":
        return ru.ru_maxrss / (1024 * 1024)
    else:
        return ru.ru_maxrss / 1024


class TimingCapture(logging.Handler):
    """
    Logging handler that captures Finished messages from timer_block to
    extract per-stage wall times without modifying production code.

    timer_block emits: "Finished <desc> (time <float>s)"
    TimerTotal emits:  "Total time for <desc>: <float>s"
    """
    FINISHED_RE = re.compile(r"Finished (.+) \(time ([\d.]+)s\)")
    TOTAL_RE = re.compile(r"Total time for (.+): ([\d.]+)s")

    def __init__(self):
        super().__init__()
        self.stages = {}

    def emit(self, record):
        msg = record.getMessage()
        m = self.FINISHED_RE.search(msg)
        if m:
            self.stages[m.group(1)] = float(m.group(2))
            return
        m = self.TOTAL_RE.search(msg)
        if m:
            self.stages[m.group(1)] = float(m.group(2))


# ---------------------------------------------------------------------------
# Core: run imputation with per-stage timing
# ---------------------------------------------------------------------------
def run_impute_timed(out_vcf_path):
    """
    Run the full Impute() pipeline, returning:
      - wall_total: total wall time
      - stages: dict of internal stage timings (from timer_block logging)
      - rss_after: peak RSS after run
      - rss_before: peak RSS before run
    """
    # Set up logging capture
    capture = TimingCapture()
    root_logger = logging.getLogger("threads_arg")
    root_logger.setLevel(logging.INFO)
    root_logger.addHandler(capture)

    rss_before = peak_rss_mb()

    t0 = time.perf_counter()
    from threads_arg.impute import Impute
    Impute(
        PANEL_VCF,
        TARGET_VCF,
        GMAP,
        MUT,
        DEMO,
        out_vcf_path,
        REGION,
    )
    wall_total = time.perf_counter() - t0

    rss_after = peak_rss_mb()

    # Clean up handler
    root_logger.removeHandler(capture)

    # Build timings dict: internal stages first, then total
    timings = {}
    for name, t in capture.stages.items():
        timings[name] = t
    timings["wall_total"] = wall_total

    return timings, rss_after, rss_before


# ---------------------------------------------------------------------------
# VCF comparison (bit-identical, ignoring date line)
# ---------------------------------------------------------------------------
def compare_vcf(expected_path, generated_path):
    """
    Compare two VCF files line-by-line. Ignores line 3 (##fileDate).
    Returns (match: bool, first_diff_line: int|None, detail: str).
    """
    with open(expected_path) as ef, open(generated_path) as gf:
        exp_lines = ef.readlines()
        gen_lines = gf.readlines()

    if len(exp_lines) != len(gen_lines):
        return False, None, f"line count mismatch: expected {len(exp_lines)}, got {len(gen_lines)}"

    for i, (exp_line, gen_line) in enumerate(zip(exp_lines, gen_lines), start=1):
        # Line 3 is the date header — skip it
        if i == 3:
            if not exp_line.startswith("##fileDate") or not gen_line.startswith("##fileDate"):
                return False, i, f"line 3 not a date header"
            continue
        if exp_line != gen_line:
            return False, i, f"expected: {exp_line[:80].rstrip()}...\n     got: {gen_line[:80].rstrip()}..."

    return True, None, "bit-identical"


# ---------------------------------------------------------------------------
# Save / load results
# ---------------------------------------------------------------------------
def save_results(tag, timings, rss_peak, rss_before, vcf_path):
    RESULTS_DIR.mkdir(parents=True, exist_ok=True)
    result = {
        "tag": tag,
        "timings": timings,
        "rss_peak_mb": rss_peak,
        "rss_before_mb": rss_before,
        "rss_delta_mb": rss_peak - rss_before,
        "vcf_path": str(vcf_path),
        "timestamp": time.strftime("%Y-%m-%d %H:%M:%S"),
    }
    out_path = RESULTS_DIR / f"{tag}.json"
    with open(out_path, "w") as f:
        json.dump(result, f, indent=2)
    print(f"  Results saved to {out_path}")
    return result


def load_results(tag):
    path = RESULTS_DIR / f"{tag}.json"
    if not path.exists():
        return None
    with open(path) as f:
        return json.load(f)


# ---------------------------------------------------------------------------
# Print helpers
# ---------------------------------------------------------------------------
def print_timings(label, timings):
    print(f"\n  {label}:")
    for stage, t in timings.items():
        print(f"    {stage:35s}  {t:8.3f}s")


def print_comparison(baseline, current):
    print("\n" + "=" * 72)
    print("  COMPARISON: baseline vs current")
    print("=" * 72)

    # Timings — show all stages from either run
    bt = baseline["timings"]
    ct = current["timings"]
    all_stages = list(dict.fromkeys(list(bt.keys()) + list(ct.keys())))

    print(f"\n  {'Stage':35s}  {'Baseline':>9s}  {'Current':>9s}  {'Speedup':>8s}")
    print(f"  {'-'*35}  {'-'*9}  {'-'*9}  {'-'*8}")
    for stage in all_stages:
        b = bt.get(stage, 0)
        c = ct.get(stage, 0)
        if b > 0 and c > 0:
            ratio = b / c
            marker = " <--" if ratio > 1.05 else (" SLOW" if ratio < 0.95 else "")
            print(f"  {stage:35s}  {b:8.3f}s  {c:8.3f}s  {ratio:7.2f}x{marker}")
        elif c > 0:
            print(f"  {stage:35s}  {'n/a':>9s}  {c:8.3f}s")
        else:
            print(f"  {stage:35s}  {b:8.3f}s  {'n/a':>9s}")

    # Memory
    print(f"\n  {'Memory':35s}  {'Baseline':>9s}  {'Current':>9s}  {'Change':>8s}")
    print(f"  {'-'*35}  {'-'*9}  {'-'*9}  {'-'*8}")
    bd = baseline["rss_delta_mb"]
    cd = current["rss_delta_mb"]
    reduction = (1 - cd / bd) * 100 if bd > 0 else 0
    sign = "-" if reduction > 0 else "+"
    print(f"  {'RSS delta':35s}  {bd:7.1f}MB  {cd:7.1f}MB  {sign}{abs(reduction):.1f}%")
    bp = baseline["rss_peak_mb"]
    cp = current["rss_peak_mb"]
    print(f"  {'RSS peak':35s}  {bp:7.1f}MB  {cp:7.1f}MB")
    print()


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------
def main():
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--compare", action="store_true",
                        help="Run as 'current' and compare against saved baseline")
    parser.add_argument("--tag", default=None,
                        help="Label for this run (default: 'baseline' or 'current')")
    parser.add_argument("--repeat", type=int, default=1,
                        help="Number of repetitions (reports best wall_total)")
    args = parser.parse_args()

    tag = args.tag or ("current" if args.compare else "baseline")

    # Validate inputs exist
    for path, label in [(PANEL_VCF, "panel"), (TARGET_VCF, "target"),
                         (GMAP, "gmap"), (MUT, "mut"), (DEMO, "demo")]:
        if not path.exists():
            print(f"ERROR: {label} not found at {path}")
            sys.exit(1)

    print(f"{'=' * 72}")
    print(f"  Imputation Benchmark — tag: {tag}")
    print(f"  Region: {REGION}  |  Repeats: {args.repeat}")
    print(f"{'=' * 72}")

    best_timings = None
    best_total = float("inf")
    rss_peak = 0
    rss_before = 0

    for rep in range(args.repeat):
        with tempfile.TemporaryDirectory() as tmpdir:
            out_vcf = Path(tmpdir) / "imputed.vcf"
            timings, rss_p, rss_b = run_impute_timed(str(out_vcf))
            wt = timings["wall_total"]

            if wt < best_total:
                best_timings = timings
                best_total = wt
                rss_peak = rss_p
                rss_before = rss_b

                # Persist the VCF from the best run
                RESULTS_DIR.mkdir(parents=True, exist_ok=True)
                vcf_out_path = RESULTS_DIR / f"{tag}.vcf"
                with open(out_vcf) as src, open(vcf_out_path, "w") as dst:
                    dst.write(src.read())

            if args.repeat > 1:
                print(f"  rep {rep+1}/{args.repeat}: {wt:.3f}s  (best: {best_total:.3f}s)")

    vcf_out_path = RESULTS_DIR / f"{tag}.vcf"

    # Print results
    print_timings(tag, best_timings)
    delta = rss_peak - rss_before
    print(f"\n  Peak RSS: {rss_peak:.1f} MB  (delta from start: {delta:.1f} MB)")

    # Bit-identity check against expected snapshot
    print(f"\n  Bit-identity vs {EXPECTED_VCF.name}...")
    match, diff_line, detail = compare_vcf(EXPECTED_VCF, vcf_out_path)
    if match:
        print(f"  PASS: bit-identical to expected snapshot")
    else:
        print(f"  FAIL: differs at line {diff_line}")
        print(f"        {detail}")

    # Save
    result = save_results(tag, best_timings, rss_peak, rss_before, vcf_out_path)

    # Compare mode
    if args.compare:
        baseline = load_results("baseline")
        if baseline is None:
            print("\n  WARNING: no baseline found. Run without --compare first.")
        else:
            print_comparison(baseline, result)

            # Bit-identity between baseline and current
            baseline_vcf = Path(baseline["vcf_path"])
            if baseline_vcf.exists():
                match2, diff2, detail2 = compare_vcf(baseline_vcf, vcf_out_path)
                if match2:
                    print("  Baseline vs current VCF: BIT-IDENTICAL")
                else:
                    print(f"  Baseline vs current VCF: DIFFERS at line {diff2}")
                    print(f"    {detail2}")
            print()

    # Summary exit
    status = "PASS" if match else "FAIL"
    print(f"  [{status}]  {tag}  {best_total:.3f}s  {rss_peak:.0f}MB")
    print()


if __name__ == "__main__":
    main()
