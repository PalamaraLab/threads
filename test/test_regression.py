import subprocess
import h5py
import numpy as np

from pathlib import Path

BASE_DIR = Path(__file__).parent.parent


def run_threads_infer():
    command = [
        "threads", "infer",
        "--pgen", str(BASE_DIR / 'example' / 'example_data.pgen'),
        "--map_gz", str(BASE_DIR / 'example' / 'example_data.map'),
        "--demography", str(BASE_DIR / 'example' / 'Ne10000.demo'),
        "--out", str(BASE_DIR / 'example' / 'example_data.threads'),
    ]

    return subprocess.run(command, cwd=BASE_DIR, capture_output=True, text=True)


def run_threads_convert():
    command = [
        "threads", "convert",
        "--threads", str(BASE_DIR / 'example' / 'example_data.threads'),
        "--argn", str(BASE_DIR / 'example' / 'example_data.argn'),
    ]

    return subprocess.run(command, cwd=BASE_DIR, capture_output=True, text=True)


def test_threads_on_example_data():
    result_infer = run_threads_infer()
    assert result_infer.returncode == 0, f"Script did not run successfully: {result_infer.stderr}"
    generated_threads = BASE_DIR / 'example' / 'example_data.threads'
    assert generated_threads.is_file()

    result_convert = run_threads_convert()
    assert result_convert.returncode == 0, f"Script did not run successfully: {result_convert.stderr}"
    generated_argn = BASE_DIR / 'example' / 'example_data.argn'
    assert generated_argn.is_file()

    reference_threads = BASE_DIR / 'test' / 'data' / 'reference.threads'
    assert reference_threads.is_file()

    reference_argn = BASE_DIR / 'test' / 'data' / 'reference.argn'
    assert reference_argn.is_file()

    with h5py.File(generated_threads, 'r') as gen_file, h5py.File(reference_threads, 'r') as ref_file:

        # Compare attributes in the root of the files - the only attribute in this file should be datetime_created
        assert gen_file.attrs.keys() == ref_file.attrs.keys()

        # Compare the names of datasets
        assert set(gen_file.keys()) == set(ref_file.keys())

        for dataset_name in gen_file.keys():
            dset_gen = gen_file[dataset_name]
            dset_ref = ref_file[dataset_name]

            print('\n\n###')
            print(dataset_name)
            print(dset_gen.shape)

            # assert dset_gen.shape == dset_ref.shape
            # assert np.allclose(dset_gen, dset_ref)
