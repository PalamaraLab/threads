import numpy as np
import h5py
import subprocess
import tempfile

from pathlib import Path

BASE_DIR = Path(__file__).parent.parent

def _run_threads_infer(outfile: Path):
    """
    Generate .threads file from example data
    """
    dir = BASE_DIR / 'example'
    command = [
        "threads", "infer",
        "--pgen", str(dir / 'example_data.pgen'),
        "--map_gz", str(dir / 'example_data.map'),
        "--demography", str(dir / 'Ne10000.demo'),
        "--out", str(outfile),
    ]

    return subprocess.run(
        command,
        cwd=BASE_DIR,
        capture_output=True,
        text=True
    )


def test_infer_consistency():
    """
    Regression check for difference in output from threads infer

    Reference data created April 19 2024 after apparently innocuous C++ changed output.
    """
    with tempfile.TemporaryDirectory() as tmpdir:
        # Regenerate threads infer output
        generated_threads = Path(tmpdir) / "test_infer_consistency.threads"
        result_infer = _run_threads_infer(generated_threads)
        assert result_infer.returncode == 0, f"Script did not run successfully: {result_infer.stderr}"
        assert generated_threads.is_file()

        # Compare against reference.threads data generated from examples dir
        reference_threads = BASE_DIR / 'test' / 'data' / 'reference.threads'
        assert reference_threads.is_file()
        with h5py.File(generated_threads, 'r') as gen_file, h5py.File(reference_threads, 'r') as ref_file:

            # Compare attributes in the root of the files - the only attribute in this file should be datetime_created
            assert gen_file.attrs.keys() == ref_file.attrs.keys()

            # Compare the names of datasets
            assert set(gen_file.keys()) == set(ref_file.keys())

            for dataset_name in gen_file.keys():
                dset_gen = gen_file[dataset_name]
                dset_ref = ref_file[dataset_name]

                assert dset_gen.shape == dset_ref.shape
                assert np.allclose(dset_gen, dset_ref)
