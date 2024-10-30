# This file is part of the Threads software suite.
# Copyright (C) 2024 Threads Developers.
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <http://www.gnu.org/licenses/>.

import numpy as np
import h5py
import tempfile

from pathlib import Path

from threads_arg.infer import threads_infer
from threads_arg.convert import threads_convert

BASE_DIR = Path(__file__).parent.parent


def _check_hdf_files_match(generated: Path, expected: Path):
    """
    Check that contents of generated hdf file matches those of expected file
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
        test_data_dir = BASE_DIR / "test" / "data"
        threads_infer(
            pgen=str(test_data_dir / "N250.pgen"),
            map_gz=str(test_data_dir / "test_data.map"),
            recombination_rate=1.3e-8,
            demography=str(test_data_dir / "Ne10000.demo"),
            mutation_rate=1.4e-8,
            query_interval=0.01,
            match_group_interval=0.5,
            mode="wgs",
            num_threads=1,
            region=None,
            max_sample_batch_size=None,
            out=str(threads_path)
        )

        # Compare against expected snapshot of threads data
        threads_expected_path = BASE_DIR / "test" / "data" / "expected_N250.threads"
        _check_hdf_files_match(threads_path, threads_expected_path)

        # Convert generated output
        argn_path = Path(tmpdir) / "test_data_snapshot_regression.argn"
        threads_convert(
            threads=str(threads_path),
            argn=str(argn_path),
            tsz=None,
            max_n=None,
            random_seed=1234,
            verify=False
        )

        # Compare against expected snapshot of argn data
        convert_expected_path = BASE_DIR / "test" / "data" / "expected_N250.argn"
        _check_hdf_files_match(argn_path, convert_expected_path)
