# This file is part of the Threads software suite.
# Copyright (C) 2024-2025 Threads Developers.
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
import pgenlib
import pytest 

from threads_arg.serialization import load_instructions

from snapshot_runners import (
    TEST_DATA_DIR
)

def _col_normalize(x):
    z = x.copy()
    mu = z.mean(axis=0, keepdims=True)
    std = z.std(axis=0, keepdims=True)
    return (z - mu) / std

def test_left_multiply():
    # Read ground truth genotypes
    pgen_path = str(TEST_DATA_DIR / "panel.pgen")
    reader = pgenlib.PgenReader(str(pgen_path).encode())
    expected_num_variants = reader.get_variant_ct()
    num_samples = reader.get_raw_sample_ct()
    expected_gt = np.empty((expected_num_variants, 2 * num_samples), dtype=np.int32)
    reader.read_alleles_range(0, expected_num_variants, expected_gt)
    gt_matrix = expected_gt.transpose()
    gt_matrix_dip = gt_matrix[::2] + gt_matrix[1::2]
    gt_matrix_norm = _col_normalize(gt_matrix)
    gt_matrix_dip_norm = _col_normalize(gt_matrix_dip)

    # Read threading instructions
    threads_path = str(TEST_DATA_DIR / "expected_infer_snapshot.threads")
    instructions = load_instructions(threads_path)

    # Random vector to multiply with
    rng = np.random.default_rng(130222)
    x_hap = rng.normal(0, 1, 2 * num_samples)
    x_dip = rng.normal(0, 1, num_samples)

    # Make sure length checks are performed
    with pytest.raises(RuntimeError):
        instructions.left_multiply(x_dip)
    with pytest.raises(RuntimeError):
        instructions.left_multiply(x_hap, diploid=True)

    # Do normal left-multiplication
    expected = x_hap @ gt_matrix
    expected_norm = x_hap @ gt_matrix_norm
    expected_dip = x_dip @ gt_matrix_dip
    expected_dip_norm = x_dip @ gt_matrix_dip_norm

    # Do threads left-multiplication and confirm results are correct
    found = instructions.left_multiply(x_hap)
    assert np.allclose(expected, found)
    found_norm = instructions.left_multiply(x_hap, normalize=True)
    assert np.allclose(expected_norm, found_norm)
    found_dip = instructions.left_multiply(x_dip, diploid=True)
    assert np.allclose(expected_dip, found_dip)
    found_dip_norm = instructions.left_multiply(x_dip, normalize=True, diploid=True)
    assert np.allclose(expected_dip_norm, found_dip_norm)
