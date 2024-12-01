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
import pandas as pd
import warnings
import logging
import pgenlib
import time

from contextlib import contextmanager

logger = logging.getLogger(__name__)

def read_map_gz(map_gz):
    """
    Reading in map file (columns 0: chrom, 1: SNP, 2: cM-pos, 3: bp)
    """
    if (map_gz[:-3] == ".gz") :
        maps = pd.read_table(map_gz, header=None, compression='gzip')
    else:
        maps = pd.read_table(map_gz, header=None)
    cm_pos = maps[2].values.astype(np.float64)
    phys_pos = maps[3].values.astype(np.float64)
    for i in range(1, len(cm_pos)):
        if cm_pos[i] <= cm_pos[i-1]:
            cm_pos[i] = cm_pos[i-1] + 1e-5
    return cm_pos, phys_pos


def interpolate_map(map_gz, pgen):
    """
    Reading in map file (format has columns [chrom, SNP, cM-pos, bp])
    """
    if (map_gz[:-3] == ".gz") :
        maps = pd.read_table(map_gz, header=None, compression='gzip', sep="\\s+")
    else:
        maps = pd.read_table(map_gz, header=None, sep="\\s+")
    cm_pos_map = maps[2].values.astype(np.float64)
    phys_pos_map = maps[3].values.astype(np.float64)
    pvar = pgen.replace("pgen", "pvar")
    bim = pgen.replace("pgen", "bim")

    physical_positions = None
    if os.path.isfile(bim):
        physical_positions = np.array(pd.read_table(bim, sep="\\s+", header=None, comment='#')[3]).astype(np.float64)
    elif os.path.isfile(pvar):
        physical_positions = np.array(pd.read_table(pvar, sep="\\s+", header=None, comment='#')[1]).astype(np.float64)
    else:
        raise RuntimeError(f"Can't find {bim} or {pvar}")

    cm_out = np.interp(physical_positions, phys_pos_map, cm_pos_map)

    if physical_positions.max() > phys_pos_map.max() or physical_positions.min() < phys_pos_map.min():
        warnings.warn("Warning: Found variants outside map range. Consider trimming input genotypes.")

    # We may get complaints in the model where the recombination rate is 0
    for i in range(1, len(cm_out)):
        if cm_out[i] <= cm_out[i-1]:
            cm_out[i] = cm_out[i-1] + 1e-5
    return cm_out, physical_positions

def read_positions_and_ids(pgen):
    pvar = pgen.replace("pgen", "pvar")
    bim = pgen.replace("pgen", "bim")

    ids, positions = [], []
    if os.path.isfile(bim):
        with open(bim, "r") as bimfile:
            for line in bimfile:
                data = line.split()
                ids.append(data[1])
                positions.append(int(data[3]))
    elif os.path.isfile(pvar):
        with open(pvar, "r") as pvarfile:
            for line in pvarfile:
                if line.startswith("#"):
                    continue
                else:
                    data = line.split()
                    ids.append(data[2])
                    positions.append(int(data[1]))
    else:
        raise RuntimeError(f"Can't find {bim} or {pvar}")
    return positions, ids

def get_map_from_bim(pgen, rho):
    pvar = pgen.replace("pgen", "pvar")
    bim = pgen.replace("pgen", "bim")
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

def iterate_pgen(pgen, callback, start_idx=None, end_idx=None, **kwargs):
    """
    Wrapper to iterate over each site in a .pgen with a callback,
    batching to reduce memory usage and read time
    """
    # Initialize read batching
    reader = pgenlib.PgenReader(pgen.encode())
    num_samples = reader.get_raw_sample_ct()
    num_sites = reader.get_variant_ct()
    if start_idx is None:
        start_idx = 0
    if end_idx is None:
        end_idx = num_sites
    M = end_idx - start_idx

    BATCH_SIZE = int(4e7 // num_samples)
    n_batches = int(np.ceil(M / BATCH_SIZE))
    # Get singleton filter for the matching step
    alleles_out = None
    phased_out = None
    i = 0
    for b in range(n_batches):
        b_start = b * BATCH_SIZE + start_idx
        b_end = min(end_idx, (b+1) * BATCH_SIZE)
        g_size = b_end - b_start
        alleles_out = np.empty((g_size, 2 * num_samples), dtype=np.int32)
        phasepresent_out = np.empty((g_size, num_samples), dtype=np.uint8)
        reader.read_alleles_and_phasepresent_range(b_start, b_end, alleles_out, phasepresent_out)
        if np.any(phasepresent_out == 0):
            raise RuntimeError("Unphased variants are currently not supported.")
        for g in alleles_out:
            callback(i, g, **kwargs)
            i += 1
    # Make sure we processed as many things as wanted to
    assert i == M
