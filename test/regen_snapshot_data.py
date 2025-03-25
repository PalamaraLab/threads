"""
Regenerate files for 'snapshot' tests.

The snapshot tests catch any regressions from a full run, they are built with
this script and later run through tests on CI. Syntax:

    python regen_snapshot_data.py <out_dir>

This is intended to be run on CI so that snapshot tests are comparing data
generated on the same architecture. After triggering the tests, the snapshot
files may be downloaded as artifacts and re-committed.
"""
import os
import sys

from pathlib import Path

from snapshot_runners import (
    run_infer_snapshot,
    run_convert_snapshot,
    run_map_snapshot,
    run_impute_snapshot
)

def regen_snapshot_files(out_dir: str):
    """
    Creates out_dir and re-runs all tests that require snapshot data.
    """
    os.mkdir(out_dir)
    out_path = Path(out_dir)

    threads_path = out_path / "expected_infer_snapshot.threads"
    run_infer_snapshot(threads_path, fit_to_data=False)

    argn_path = out_path / "expected_convert_snapshot.argn"
    run_convert_snapshot(threads_path, argn_path)

    mut_path = out_path / "expected_mapping_snapshot.mut"
    run_map_snapshot(argn_path, mut_path)

    vcf_path = out_path / "expected_impute_snapshot.vcf"
    run_impute_snapshot(mut_path, vcf_path)

    # Re-run infer and convert with fit_to_data
    threads_ftd_path = out_path / "expected_infer_fit_to_data_snapshot.threads"
    run_infer_snapshot(threads_ftd_path, fit_to_data=True)

    argn_ftd_path = out_path / "expected_convert_fit_to_data_snapshot.argn"
    run_convert_snapshot(threads_ftd_path, argn_ftd_path)


if __name__ == "__main__":
    if len(sys.argv) < 2:
        raise RuntimeError("Specify output dir as first argument")

    regen_snapshot_files(sys.argv[1])
