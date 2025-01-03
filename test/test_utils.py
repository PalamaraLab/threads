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

import pytest

from threads_arg.utils import parse_region_string

def test_parse_region_string():
    # Check valid cases
    assert parse_region_string("1:2-3") == (1, 2, 3)
    assert parse_region_string("chr1:2-3") == (1, 2, 3)
    assert parse_region_string("2-3") == (None, 2, 3)
    assert parse_region_string("22:2000-3000") == (22, 2000, 3000)
    assert parse_region_string("chr22:2000-3000") == (22, 2000, 3000)
    assert parse_region_string("2000-3000") == (None, 2000, 3000)

    # Check invalid cases
    with pytest.raises(Exception) as exc_info:
        parse_region_string("totallywrong")
    assert exc_info.value.args[0] == "Invalid region string 'totallywrong'"

    with pytest.raises(Exception) as exc_info:
        parse_region_string("fish1:2-3")
    assert exc_info.value.args[0] == "Invalid region string 'fish1:2-3'"

    with pytest.raises(Exception) as exc_info:
        parse_region_string("chr0:2-3")
    assert exc_info.value.args[0] == "Invalid chromosome 0 not between 1 and 22"

    with pytest.raises(Exception) as exc_info:
        parse_region_string("chr23:2-3")
    assert exc_info.value.args[0] == "Invalid chromosome 23 not between 1 and 22"

    with pytest.raises(Exception) as exc_info:
        parse_region_string("chr20:7000-5000")
    assert exc_info.value.args[0] == "Invalid range: end 5000 less than start 7000"
