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

def decompress_threads(threads):
    """
    Read in threading instructions for the purpose of converting to tskit/arg-needle format.
    """
    f = h5py.File(threads, "r")

    samples, thread_starts, mut_starts = f["samples"][:, 0], f["samples"][:, 1], f["samples"][:, 2]
    positions = f['positions'][...]
    flat_ids, flat_bps = f['thread_targets'][:, :-1], f['thread_targets'][:, -1]
    flat_ages = f['thread_ages'][...]
    flat_hets = f['het_sites'][...]
    arg_range = f['arg_range'][...]

    threading_instructions = []
    mutations = []
    for i, (start, start_mut) in enumerate(zip(thread_starts, mut_starts)):
        if i == len(thread_starts) - 1:
            ids = flat_ids[start:]
            bps = flat_bps[start:]
            ages = flat_ages[start:]
            muts = flat_hets[start_mut:]
        else:
            ids = flat_ids[start:thread_starts[i + 1]]
            bps = flat_bps[start:thread_starts[i + 1]]
            ages = flat_ages[start:thread_starts[i + 1]]
            muts = flat_hets[start_mut:mut_starts[i + 1]]
        threading_instructions.append((bps, ids, ages))
        mutations.append(muts)
    
    return {
        "threads": threading_instructions,
        "samples": samples,
        "positions": positions,
        "arg_range": arg_range,
        "mutations": mutations
    }


def slice_threading_instructions(threading_data, start_pos, end_pos):
    """
    Slice threading instructions to only include data within the region
    start_pos-end_pos (both inclusive)
    TODO: make threading instructions be a (ideally c++-side) class object and not a random dict.
    """
    if start_pos >= end_pos:
        raise ValueError(f"Need start_pos < end_pos, found start_pos={start_pos} and end_pos={end_pos}")


    old_range = threading_data["arg_range"]
    if start_pos < old_range[0] or start_pos >= old_range[1]:
        raise ValueError(f"Need start_pos within arg range, found start_pos={start_pos} and range {old_range[0]}-{old_range[1]}")
    if end_pos <= old_range[0] or end_pos > old_range[1]:
        raise ValueError(f"Need end_pos within arg range, found end_pos={end_pos} and range {old_range[0]}-{old_range[1]}")

    positions = threading_data["positions"]

    new_threading_instructions = []
    new_mutations = []
    new_arg_range = [start_pos, end_pos]
    new_samples = threading_data["samples"]
    new_positions = positions[(start_pos <= positions) & (positions <= end_pos)]

    for muts, (bps, ids, ages) in zip(threading_data["mutations"], threading_data["threads"]):
        new_mutations.append(muts[(start_pos <= muts) & (muts <= end_pos)])

        if len(bps) == 0:
            new_threading_instructions.append((np.array([]), np.array([]), np.array([])))
            continue

        first_idx = np.nonzero(bps <= start_pos)[0][-1]
        if bps[-1] <= end_pos:
            last_idx = len(bps)
        else:
            last_idx = np.nonzero(end_pos <= bps)[0][0]
        new_bps = np.concatenate([np.array([start_pos]), bps[(first_idx + 1):last_idx]])
        new_ids = ids[first_idx:last_idx]
        new_ages = ages[first_idx:last_idx]
        new_threading_instructions.append((new_bps, new_ids, new_ages))

    return {
        "threads": new_threading_instructions,
        "samples": new_samples,
        "positions": new_positions,
        "arg_range": new_arg_range,
        "mutations": new_mutations
    }


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
