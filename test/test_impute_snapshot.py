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

import tempfile
import time

from pathlib import Path

from threads_arg.impute import Impute

BASE_DIR = Path(__file__).parent.parent


def _check_vcf_files_match(expected_vcf_path, generated_vcf_path):
    expected = open(expected_vcf_path)
    generated = open(generated_vcf_path)

    for i, lhs in enumerate(expected):
        line = i + 1
        rhs = generated.readline()
        if line == 3:
            assert(lhs[0:10] == "##fileDate")
        else:
            assert(lhs == rhs)


def test_impute_snapshot_regression():
    """
    Regression check for difference in output from threads impute
    """
    with tempfile.TemporaryDirectory() as tmp_dir:
        generated_vcf_path = Path(tmp_dir) / "test_imputed.vcf"
        data_dir = BASE_DIR / "test" / "data"

        start_time = time.perf_counter()
        Impute(
            data_dir / "panel.vcf.gz",
            data_dir / "target.vcf.gz",
            data_dir / "gmap_04.map",
            data_dir / "mapping.mut",
            data_dir / "CEU_unscaled.demo",
            generated_vcf_path,
            "1:400000-600000"
        )

        end_time = time.perf_counter()
        print(f"Impute ran in {end_time - start_time}s")

        _check_vcf_files_match(data_dir / "imputed.vcf", generated_vcf_path)
