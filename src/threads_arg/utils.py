import os
import numpy as np
import h5py
import pandas as pd
import warnings

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
