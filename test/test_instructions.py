# This file is part of the Threads software suite.
# Copyright (C) 2025 Threads Developers.
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

import pickle

from threads_arg import ThreadingInstructions

def test_threading_instructions_pickle():
    ti = ThreadingInstructions(
        [[21, 22, 23]], # Not valid data; just used for comparison below
        [[31, 32, 33]],
        [[41, 42, 43]],
        [[0, 1, 2, 3, 4]],
        [10, 20, 30, 40, 50],
        10, 50
    )

    ti_pickle = pickle.dumps(ti)
    ti_restored = pickle.loads(ti_pickle)

    assert ti_restored.all_starts() == ti.all_starts()
    assert ti_restored.all_tmrcas() == ti.all_tmrcas()
    assert ti_restored.all_targets() == ti.all_targets()
    assert ti_restored.all_mismatches() == ti.all_mismatches()
    assert ti_restored.positions == ti.positions
    assert ti_restored.num_sites == ti.num_sites
    assert ti_restored.num_samples == ti.num_samples
    assert ti_restored.start == ti.start
    assert ti_restored.end == ti.end
