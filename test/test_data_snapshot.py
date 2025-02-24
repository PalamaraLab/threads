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
import pgenlib
import arg_needle_lib

from pathlib import Path

from threads_arg.infer import threads_infer
from threads_arg.convert import threads_convert
from threads_arg.serialization import load_instructions
from threads_arg import GenotypeIterator

BASE_DIR = Path(__file__).parent.parent


def assert_shape_msg(generated, expected, dataset_name, gen_item, exp_item):
    msg = f"Difference in dataset shape '{dataset_name}' dataset\n"
    msg += f"Generated {gen_item.shape} ({generated}):\n{np.array(gen_item)}\n"
    msg += f"Expected {exp_item.shape} ({expected}):\n{np.array(exp_item)}\n"
    return msg


def assert_allclose_msg(generated, expected, dataset_name, gen_item, exp_item):
    msg = f"Detected diff in '{dataset_name}' dataset\n"
    msg += f"Generated ({generated}):\n{np.array(gen_item)}\n"
    msg += f"Expected ({expected}):\n{np.array(exp_item)}\n"
    return msg


def _check_hdf_items_match(generated_file, expected_file, generated_items, expected_items, path):
    """
    Check that datasets or groups match. Recurses on groups.
    """
    for name in generated_items.keys():
        gen_item = generated_items[name]
        exp_item = expected_items[name]
        item_path = f"{path}/{name}"
        if isinstance(gen_item, h5py.Group):
            # Assert that sibling group exists and recurse into child
            assert isinstance(exp_item, h5py.Group)
            _check_hdf_items_match(generated_file, expected_file, gen_item, exp_item, item_path)
        else:
            # Assert that shape same and values close enough
            assert gen_item.shape == exp_item.shape, assert_shape_msg(generated_file, expected_file, name, gen_item, exp_item)
            assert np.allclose(gen_item, exp_item), assert_allclose_msg(generated_file, expected_file, name, gen_item, exp_item)


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

        _check_hdf_items_match(generated, expected, gen_file, exp_file, f"{expected.name}: ")


def _check_compression_is_correct(threads: Path, pgen: Path):
    """
    Check that the threading instructions match the input pgen
    """
    reader = pgenlib.PgenReader(str(pgen).encode())
    expected_num_variants = reader.get_variant_ct()
    num_samples = reader.get_raw_sample_ct()

    expected_gt = np.empty((expected_num_variants, 2 * num_samples), dtype=np.int32)
    reader.read_alleles_range(0, expected_num_variants, expected_gt)

    instructions = load_instructions(threads)
    gt_it = GenotypeIterator(instructions)
    found_gt = np.empty((expected_num_variants, 2 * num_samples), dtype=np.int32)

    i = 0
    while gt_it.has_next_genotype():
        found_gt[i] = np.array(gt_it.next_genotype())
        i += 1

    assert (found_gt == expected_gt).all()


def _check_argn_mutations_fit_to_data(argn: Path, pgen: Path):
    reader = pgenlib.PgenReader(str(pgen).encode())
    expected_num_variants = reader.get_variant_ct()
    num_samples = reader.get_raw_sample_ct()

    expected_gt = np.empty((expected_num_variants, 2 * num_samples), dtype=np.int32)
    reader.read_alleles_range(0, expected_num_variants, expected_gt)

    arg = arg_needle_lib.deserialize_arg(str(argn))
    arg.populate_children_and_roots()
    assert len(arg.mutations()) == expected_num_variants

    argn_gt = arg_needle_lib.get_mutations_matrix(arg)
    assert (argn_gt == expected_gt).all()


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
            pgen=str(test_data_dir / "panel.pgen"),
            map=str(test_data_dir / "gmap_02.map"),
            recombination_rate=1.3e-8,
            demography=str(test_data_dir / "CEU_unscaled.demo"),
            mutation_rate=1.4e-8,
            query_interval=0.01,
            match_group_interval=0.5,
            mode="wgs",
            num_threads=1,
            region=None,
            fit_to_data=False,
            allele_ages=None,
            max_sample_batch_size=None,
            out=str(threads_path)
        )

        # Compare genotypes are correctly compressed
        _check_compression_is_correct(threads_path, str(test_data_dir / "panel.pgen"))

        # Compare against expected snapshot of threads data
        threads_expected_path = test_data_dir / "arg.threads"
        _check_hdf_files_match(threads_path, threads_expected_path)

        # Convert generated output
        argn_path = Path(tmpdir) / "test_data_snapshot_regression.argn"
        threads_convert(
            threads=str(threads_path),
            argn=str(argn_path),
            tsz=None
        )

        # Compare against expected snapshot of argn data
        convert_expected_path = test_data_dir / "arg.argn"
        _check_hdf_files_match(argn_path, convert_expected_path)


def test_fit_to_data_snapshot_regression():
    """
    Regression check for difference in output from threads infer and convert using
    the --fit-to-data flag
    """
    with tempfile.TemporaryDirectory() as tmpdir:
        # Regenerate threads infer output
        threads_path = Path(tmpdir) / "test_fit_to_data_snapshot_regression.threads"
        test_data_dir = BASE_DIR / "test" / "data"
        threads_infer(
            pgen=str(test_data_dir / "panel.pgen"),
            map=str(test_data_dir / "gmap_02.map"),
            recombination_rate=1.3e-8,
            demography=str(test_data_dir / "CEU_unscaled.demo"),
            mutation_rate=1.4e-8,
            query_interval=0.01,
            match_group_interval=0.5,
            mode="wgs",
            num_threads=1,
            region=None,
            fit_to_data=True,
            allele_ages=None,
            max_sample_batch_size=None,
            out=str(threads_path)
        )

        # Compare genotypes are correctly compressed
        _check_compression_is_correct(threads_path, str(test_data_dir / "panel.pgen"))

        # Compare against expected snapshot of threads data
        threads_expected_path = test_data_dir / "arg_consistent.threads"
        _check_hdf_files_match(threads_path, threads_expected_path)

        # Convert generated output
        argn_path = Path(tmpdir) / "test_fit_to_data_snapshot_regression.argn"
        threads_convert(
            threads=str(threads_path),
            argn=str(argn_path),
            tsz=None,
            add_mutations=True
        )

        # Compare against expected snapshot of argn data
        convert_expected_path = test_data_dir / "arg_consistent.argn"
        _check_hdf_files_match(argn_path, convert_expected_path)

        # Put this back in here once arg_needle_lib-dev #19 and #20 have been resolved
        # _check_argn_mutations_fit_to_data(argn_path, str(test_data_dir / "panel.pgen"))

