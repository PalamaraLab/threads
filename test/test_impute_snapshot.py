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

import tempfile
import time

from pathlib import Path

from snapshot_runners import (
    run_map_snapshot,
    run_impute_snapshot,
    TEST_DATA_DIR
)

BASE_DIR = Path(__file__).parent.parent


def _line_check_eq(_line_no, lhs, rhs):
    return lhs == rhs


def _line_check_eq_ignore_date(line_no, lhs, rhs):
    if line_no == 3:
        return lhs[0:10] == "##fileDate"
    else:
        return lhs == rhs


def _check_files_match(expected_path, generated_path, line_compare_fn):
    expected = open(expected_path)
    generated = open(generated_path)

    for i, lhs in enumerate(expected):
        line_no = i + 1
        rhs = generated.readline()
        if not line_compare_fn(line_no, lhs, rhs):
            assert lhs == rhs, f"line no {line_no} differs between expected {expected_path} and generated {generated_path}"


def test_map_snapshot_regression():
    """
    Regression check for for deterministic mapping generation
    """
    with tempfile.TemporaryDirectory() as tmpdir:
        # Pre-generated argn snapshot used as valid mapping input
        generated_mut_path = Path(tmpdir) / "mapping.mut"
        start_time = time.perf_counter()
        run_map_snapshot(
            TEST_DATA_DIR / "expected_convert_snapshot.argn",
            generated_mut_path
        )
        end_time = time.perf_counter()
        print(f"Map ran in {end_time - start_time}s")

        _check_files_match(
            expected_path=str(TEST_DATA_DIR / "expected_mapping_snapshot.mut"),
            generated_path=str(generated_mut_path),
            line_compare_fn=_line_check_eq
        )


def test_impute_snapshot_regression():
    """
    Regression check for difference in output from threads impute
    """
    with tempfile.TemporaryDirectory() as tmpdir:
        # Pre-generated mut snapshot used as valid imputation input
        generated_vcf_path = Path(tmpdir) / "imputed.vcf"
        start_time = time.perf_counter()
        run_impute_snapshot(
            TEST_DATA_DIR / "expected_mapping_snapshot.mut",
            generated_vcf_path
        )
        end_time = time.perf_counter()
        print(f"Impute ran in {end_time - start_time}s")

        # The files should be identical except one line with new date
        _check_files_match(
            expected_path=str(TEST_DATA_DIR / "expected_impute_snapshot.vcf"),
            generated_path=str(generated_vcf_path),
            line_compare_fn=_line_check_eq_ignore_date
        )
