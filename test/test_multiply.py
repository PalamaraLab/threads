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

"""
Tests for genotype-matrix multiplication.

Dense multiply is the reference baseline. RLE (sparse baseline) and DAG
(fast sublinear multiply, following the arg-needle-lib approach) are
validated against it.
"""

import numpy as np
import pgenlib
import pytest

from threads_arg.serialization import load_instructions
from snapshot_runners import TEST_DATA_DIR


# ---------------------------------------------------------------------------
# Fixtures
# ---------------------------------------------------------------------------

@pytest.fixture(scope="module")
def ground_truth():
    """Load ground-truth genotype matrix from pgen."""
    pgen_path = str(TEST_DATA_DIR / "panel.pgen")
    reader = pgenlib.PgenReader(pgen_path.encode())
    n_var = reader.get_variant_ct()
    n_samp = reader.get_raw_sample_ct()
    gt = np.empty((n_var, 2 * n_samp), dtype=np.int32)
    reader.read_alleles_range(0, n_var, gt)
    return gt.T  # (2*n_samp, n_var)


@pytest.fixture(scope="module")
def instructions():
    """Load threading instructions and prepare RLE + DAG."""
    path = str(TEST_DATA_DIR / "expected_infer_snapshot.threads")
    ti = load_instructions(path)
    ti.prepare_rle_multiply()
    ti.prepare_dag_multiply(num_chunks=1)
    return ti


def _col_normalize(x):
    mu = x.mean(axis=0, keepdims=True)
    std = x.std(axis=0, keepdims=True)
    return (x - mu) / std


# ---------------------------------------------------------------------------
# Parametrized tests: each method vs dense baseline
# ---------------------------------------------------------------------------

METHODS = ["dense", "rle", "dag"]
MODES = [
    pytest.param(False, False, id="haploid"),
    pytest.param(True, False, id="diploid"),
    pytest.param(False, True, id="normalize"),
    pytest.param(True, True, id="diploid+normalize"),
]


def _right_multiply(ti, method, x, diploid, normalize):
    if method == "dense":
        return ti.right_multiply(x, diploid=diploid, normalize=normalize)
    elif method == "rle":
        return ti.right_multiply_rle(x, diploid=diploid, normalize=normalize)
    elif method == "dag":
        return ti.right_multiply_dag(x, diploid=diploid, normalize=normalize)


def _left_multiply(ti, method, x, diploid, normalize):
    if method == "dense":
        return ti.left_multiply(x, diploid=diploid, normalize=normalize)
    elif method == "rle":
        return ti.left_multiply_rle(x, diploid=diploid, normalize=normalize)
    elif method == "dag":
        return ti.left_multiply_dag(x, diploid=diploid, normalize=normalize)


@pytest.mark.parametrize("method", METHODS)
@pytest.mark.parametrize("diploid,normalize", MODES)
def test_right_multiply(ground_truth, instructions, method, diploid, normalize):
    """G @ x must match dense numpy baseline."""
    gt = ground_truth  # (n_hap, m)
    m = gt.shape[1]
    rng = np.random.default_rng(42)
    x = rng.normal(size=m)

    if diploid:
        ref = _col_normalize(gt[::2] + gt[1::2]) @ x if normalize else (gt[::2] + gt[1::2]) @ x
    else:
        ref = _col_normalize(gt.astype(float)) @ x if normalize else gt @ x

    result = _right_multiply(instructions, method, x, diploid, normalize)
    assert np.allclose(ref, result, atol=1e-10), f"{method} right {diploid=} {normalize=}"


@pytest.mark.parametrize("method", METHODS)
@pytest.mark.parametrize("diploid,normalize", MODES)
def test_left_multiply(ground_truth, instructions, method, diploid, normalize):
    """G.T @ x must match dense numpy baseline."""
    gt = ground_truth  # (n_hap, m)
    n_hap = gt.shape[0]
    n_dip = n_hap // 2
    rng = np.random.default_rng(99)

    if diploid:
        x = rng.normal(size=n_dip)
        G = gt[::2] + gt[1::2]
        ref = _col_normalize(G.astype(float)).T @ x if normalize else x @ G
    else:
        x = rng.normal(size=n_hap)
        ref = _col_normalize(gt.astype(float)).T @ x if normalize else x @ gt

    result = _left_multiply(instructions, method, x, diploid, normalize)
    assert np.allclose(ref, result, atol=1e-10), f"{method} left {diploid=} {normalize=}"


# ---------------------------------------------------------------------------
# Batch multiply (RLE + DAG)
# ---------------------------------------------------------------------------

@pytest.mark.parametrize("method", ["rle", "dag"])
def test_right_batch(ground_truth, instructions, method):
    """Batch right multiply matches column-by-column."""
    m = ground_truth.shape[1]
    k = 3
    rng = np.random.default_rng(7)
    X = rng.normal(size=(m, k))

    singles = np.column_stack([
        _right_multiply(instructions, method, X[:, j], False, False) for j in range(k)
    ])
    batch = _right_multiply(instructions, method, X, False, False)
    assert np.allclose(singles, batch, atol=1e-10)


@pytest.mark.parametrize("method", ["rle", "dag"])
def test_left_batch(ground_truth, instructions, method):
    """Batch left multiply matches column-by-column."""
    n_hap = ground_truth.shape[0]
    k = 3
    rng = np.random.default_rng(8)
    X = rng.normal(size=(n_hap, k))

    singles = np.column_stack([
        _left_multiply(instructions, method, X[:, j], False, False) for j in range(k)
    ])
    batch = _left_multiply(instructions, method, X, False, False)
    assert np.allclose(singles, batch, atol=1e-10)


# ---------------------------------------------------------------------------
# Input validation
# ---------------------------------------------------------------------------

def test_length_checks(instructions):
    """Wrong-length inputs raise RuntimeError."""
    bad = np.ones(instructions.num_sites + 1)
    with pytest.raises(RuntimeError):
        instructions.left_multiply(bad)
    with pytest.raises(RuntimeError):
        instructions.right_multiply(bad)
