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

import os
import numpy as np
import h5py
import pandas as pd
import warnings
import logging
import time

from contextlib import contextmanager

logger = logging.getLogger(__name__)


def decompress_threads(threads):
    f = h5py.File(threads, "r")

    samples, thread_starts = f["samples"][:, 0], f["samples"][:, 1]
    positions = f['positions'][...]
    flat_ids, flat_bps = f['thread_targets'][:, :-1], f['thread_targets'][:, -1]
    flat_ages = f['thread_ages'][...]
    try:
        arg_range = f['arg_range'][...]
    except KeyError:
        arg_range = [np.nan, np.nan]

    threading_instructions = []
    for i, start in enumerate(thread_starts):
        if i == len(thread_starts) - 1:
            ids = flat_ids[start:]
            bps = flat_bps[start:]
            ages = flat_ages[start:]
        else:
            ids = flat_ids[start:thread_starts[i + 1]]
            bps = flat_bps[start:thread_starts[i + 1]]
            ages = flat_ages[start:thread_starts[i + 1]]
        threading_instructions.append((bps, ids, ages))
    return {
        "threads": threading_instructions,
        "samples": samples,
        "positions": positions,
        "arg_range": arg_range
    }


def read_map_file(map_file):
    """
    Reading in map file for Li-Stephens using genetic maps in the SHAPEIT format
    """
    maps = pd.read_table(map_file, sep=r"\s+")
    cm_pos = maps.cM.values.astype(np.float64)
    phys_pos = maps.pos.values.astype(np.float64)

    # Add an epsilon value to avoid div zero errors. A few decimal places is
    # enough to distinguish each cM, so 1e-5 does not skew results.
    for i in range(1, len(cm_pos)):
        if cm_pos[i] <= cm_pos[i-1]:
            cm_pos[i] = cm_pos[i-1] + 1e-5
    return phys_pos, cm_pos


def interpolate_map(map_file, pgen_file):
    """
    Read map file in SHAPEIT format to interpolate read pgen positions to cM
    """
    phys_pos, cm_pos = read_map_file(map_file)

    pvar = pgen_file.replace("pgen", "pvar")
    bim = pgen_file.replace("pgen", "bim")
    physical_positions = None
    if os.path.isfile(bim):
        physical_positions = np.array(pd.read_table(bim, sep="\\s+", header=None, comment='#')[3]).astype(np.float64)
    elif os.path.isfile(pvar):
        physical_positions = np.array(pd.read_table(pvar, sep="\\s+", header=None, comment='#')[1]).astype(np.float64)
    else:
        raise RuntimeError(f"Can't find {bim} or {pvar}")

    cm_out = np.interp(physical_positions, phys_pos, cm_pos)

    if physical_positions.max() > phys_pos.max() or physical_positions.min() < phys_pos.min():
        warnings.warn("Warning: Found variants outside map range. Consider trimming input genotypes.")

    # We may get complaints in the model where the recombination rate is 0
    for i in range(1, len(cm_out)):
        if cm_out[i] <= cm_out[i-1]:
            cm_out[i] = cm_out[i-1] + 1e-5
    return cm_out, physical_positions


def get_map_from_bim(pgen_file, rho):
    pvar = pgen_file.replace("pgen", "pvar")
    bim = pgen_file.replace("pgen", "bim")
    cm_out = None
    physical_positions = None
    if os.path.isfile(bim):
        physical_positions = np.array(pd.read_table(bim, sep="\\s+", header=None, comment='#')[3]).astype(int)
        cm_out = rho * 100 * physical_positions
    elif os.path.isfile(pvar):
        physical_positions = np.array(pd.read_table(pvar, sep="\\s+", header=None, comment='#')[1]).astype(int)
        cm_out = rho * 100 * physical_positions
    else:
        raise RuntimeError(f"Can't find {bim} or {pvar}")

    for i in range(1, len(cm_out)):
        if cm_out[i] <= cm_out[i-1]:
            cm_out[i] = cm_out[i-1] + 1e-5
    return cm_out, physical_positions


def parse_demography(demography):
    d = pd.read_table(demography, sep="\\s+", header=None)
    return list(d[0]), list(d[1])


def split_list(list, n):
    """Yield n number of sequential chunks from l."""
    sublists = []
    d, r = divmod(len(list), n)
    for i in range(n):
        si = (d+1)*(i if i < r else r) + d*(0 if i < r else i - r)
        sublists.append(list[si:si+(d+1 if i < r else d)])
    return sublists


@contextmanager
def timer_block(desc: str, print_start: bool=True):
    """
    Context manager to log a description and time spend inside `with` block.

    By default this prints description at block start and end. Set print_start
    to false to disable the former print. This is neater for very quick blocks.

    Example usage:
        with timer_block("expensive op"):
            sleep(1)
        # Logger info happens here
    """
    if print_start:
        logger.info(f"Starting {desc}...")
    start_time = time.time()
    yield
    end_time = time.time()
    duration = end_time - start_time
    logger.info(f"Finished {desc} (time {duration:.3f}s)")


class TimerTotal:
    """
    Tally up repeated timer measurements using `with` blocks to get a string
    summary of total time.

    Example usage:
        tt = TimerTotal("foo")
        for _ in range(5):
            with tt:
                time.sleep(.4)
        print(tt)
    """
    def __init__(self, desc: str):
        self.desc = desc
        self.durations = []
        self._start_time = None

    def __enter__(self):
        self._start_time = time.time()

    def __exit__(self, _exc_type, _exc_val, _exc_tb):
        end_time = time.time()
        self.durations.append(end_time - self._start_time)
        self._start_time = None

    def __str__(self):
        total = sum(self.durations)
        return f"Total time for {self.desc}: {total:.3f}s"
