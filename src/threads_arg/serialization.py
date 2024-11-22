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

import h5py
import numpy as np
from datetime import datetime
from threads_arg import ThreadingInstructions

def serialize_instructions(instructions, out):
    num_threads = instructions.num_samples
    num_sites = instructions.num_sites
    positions = instructions.positions
    region_start = instructions.start
    region_end = instructions.end
    samples = list(range(num_threads))

    all_starts = instructions.all_starts()
    all_targets = instructions.all_targets()
    all_tmrcas = instructions.all_tmrcas()
    all_mismatches = instructions.all_mismatches()

    thread_starts = np.cumsum([0] + [len(starts) for starts in all_starts[:-1]])
    mut_starts = np.cumsum([0] + [len(mismatches) for mismatches in all_mismatches[:-1]])

    flat_starts = [start for starts in all_starts for start in starts]
    flat_tmrcas = [tmrca for tmrcas in all_tmrcas for tmrca in tmrcas]
    flat_targets = [target for targets in all_targets for target in targets]
    flat_mismatches = [mismatch for mismatches in all_mismatches for mismatch in mismatches]

    num_stitches = len(flat_starts)
    num_mutations = len(flat_mismatches)

    f = h5py.File(out, "w")
    f.attrs['datetime_created'] = datetime.now().isoformat()

    compression_opts = 9
    dset_samples = f.create_dataset("samples", (num_threads, 3), dtype=int, compression='gzip',
                                    compression_opts=compression_opts)
    dset_pos = f.create_dataset("positions", (num_sites), dtype=int, compression='gzip',
                                    compression_opts=compression_opts)
    # First L columns are random samples for imputation
    dset_targets = f.create_dataset("thread_targets", (num_stitches, 2), dtype=int, compression='gzip',
                                    compression_opts=compression_opts)
    dset_ages = f.create_dataset("thread_ages", (num_stitches), dtype=np.double, compression='gzip',
                                    compression_opts=compression_opts)
    dset_het_s = f.create_dataset("het_sites", (num_mutations), dtype=int, compression='gzip',
                                    compression_opts=compression_opts)
    dset_range = f.create_dataset("arg_range", (2), dtype=np.double, compression='gzip',
                                  compression_opts=compression_opts)

    dset_samples[:, 0] = samples
    dset_samples[:, 1] = thread_starts
    dset_samples[:, 2] = mut_starts

    dset_targets[:, 0] = flat_targets
    dset_targets[:, 1] = flat_starts

    dset_pos[:] = positions
    dset_ages[:] = flat_tmrcas
    dset_het_s[:] = flat_mismatches
    dset_range[:] = [region_start, region_end]

    f.close()

def load_instructions(threads):
    f = h5py.File(threads, "r")

    _, thread_starts, het_starts = f["samples"][:, 0], f["samples"][:, 1], f["samples"][:, 2]
    positions = f['positions'][...]
    flat_targets, flat_starts = f['thread_targets'][:, 0], f['thread_targets'][:, -1]
    flat_tmrcas = f['thread_ages'][...]
    flat_mismatches = f['het_sites'][...]

    try:
        arg_range = f['arg_range'][...]
    except KeyError:
        arg_range = [np.nan, np.nan]

    region_start = int(arg_range[0])
    region_end = int(arg_range[1])

    starts = []
    targets = []
    tmrcas = []
    mismatches = []
    for i, (start, het_start) in enumerate(zip(thread_starts, het_starts)):
        if i == len(thread_starts) - 1:
            targets.append(flat_targets[start:].tolist())
            starts.append(flat_starts[start:].tolist())
            tmrcas.append(flat_tmrcas[start:].tolist())
            mismatches.append(flat_mismatches[het_start:].tolist())
        else:
            targets.append(flat_targets[start:thread_starts[i + 1]].tolist())
            starts.append(flat_starts[start:thread_starts[i + 1]].tolist())
            tmrcas.append(flat_tmrcas[start:thread_starts[i + 1]].tolist())
            mismatches.append(flat_mismatches[het_start:het_starts[i + 1]].tolist())
    return ThreadingInstructions(starts, tmrcas, targets, mismatches, positions.astype(int).tolist(), region_start, region_end)
