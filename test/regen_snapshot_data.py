import os
import sys
import threads_arg


def regen_snapshot_files(out_dir: str):
    """
    Creates out_dir and re-runs all tests that require snapshot data.

    This is intended to be run on CI so that snapshot tests are comparing data
    generated on the same architecture. After triggering the tests, the snapshot
    files may be downloaded as artifacts and re-committed.
    """
    os.mkdir(out_dir)

    out_files = [
        "expected_infer_snapshot.argn",
        "expected_convert_snapshot.argn",
        "expected_mapping_snapshot.mut",
        "expected_imputed_snapshot.vcf"
    ]

    for out_file in out_files:
        with open(f"{out_dir}/{out_file}", "wt") as f:
            f.write(f"Dummy data for {out_file}")


if __name__ == "__main__":
    if len(sys.argv) < 2:
        raise RuntimeError("Specify output dir as first argument")

    regen_snapshot_files(sys.argv[1])
