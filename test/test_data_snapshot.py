import numpy as np
import h5py
import subprocess
import tempfile

from pathlib import Path

BASE_DIR = Path(__file__).parent.parent


def _run_threads_arg(in_pgen_filename: str, out_threads_file: Path):
    """
    Generate .threads file from test data
    in_pgen_filename is relative to test data dir
    """
    dir = BASE_DIR / "test" / "data"
    command = [
        "python3",
        "-m", "threads_arg",
        "infer",
        "--pgen", str(dir / in_pgen_filename),
        "--map_gz", str(dir / "test_data.map"),
        "--demography", str(dir / "Ne10000.demo"),
        "--out", str(out_threads_file),
    ]

    return subprocess.run(
        command,
        cwd=BASE_DIR,
        capture_output=True,
        text=True
    )


def _run_threads_convert(in_threads_file: Path, out_argn_file: Path):
    """
    Convert .threads file to argn
    """
    command = [
        "python3",
        "-m", "threads_arg",
        "convert",
        "--threads", str(in_threads_file),
        "--argn", str(out_argn_file),
        "--random-seed", "1234"
    ]

    return subprocess.run(
        command,
        cwd=BASE_DIR,
        capture_output=True,
        text=True
    )


def _check_hdf_files_match(generated: Path, expected: Path):
    """
    Check that contents of generated hdf file matche those of expected file
    """
    assert generated.is_file()
    assert expected.is_file()
    with h5py.File(generated, "r") as gen_file, h5py.File(expected, "r") as exp_file:

        # Compare attributes in the root of the files - the only attribute in this file should be datetime_created
        assert gen_file.attrs.keys() == exp_file.attrs.keys()

        # Compare the names of datasets
        assert set(gen_file.keys()) == set(exp_file.keys())

        for dataset_name in gen_file.keys():
            dset_gen = gen_file[dataset_name]
            dset_ref = exp_file[dataset_name]

            assert dset_gen.shape == dset_ref.shape
            assert np.allclose(dset_gen, dset_ref)


def test_data_snapshot_regression():
    """
    Regression check for difference in output from threads infer and convert

    This was added after an upgrade to pybind caused a change in generated data
    and is being used to determine which parts of code need further testing.
    """
    with tempfile.TemporaryDirectory() as tmpdir:
        # Regenerate threads infer output
        threads_path = Path(tmpdir) / "test_data_snapshot_regression.threads"
        infer_result = _run_threads_arg("N250.pgen", threads_path)
        assert infer_result.returncode == 0, f"threads infer did not run successfully: {infer_result.stderr}"

        # Compare against reference.threads data generated from examples dir
        threads_expected_path = BASE_DIR / "test" / "data" / "expected_N250.threads"
        _check_hdf_files_match(threads_path, threads_expected_path)

        # Convert generated output
        argn_path = Path(tmpdir) / "test_data_snapshot_regression.argn"
        convert_result = _run_threads_convert(threads_path, argn_path)
        assert convert_result.returncode == 0, f"threads convert did not run successfully: {infer_result.stderr}"
        convert_expected_path = BASE_DIR / "test" / "data" / "expected_N250.argn"
        _check_hdf_files_match(argn_path, convert_expected_path)
