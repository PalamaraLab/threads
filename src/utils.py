import numpy as np
import h5py
import pandas as pd

def decompress_threads(threads):
    f = h5py.File(threads, "r")

    samples, thread_starts = f["samples"][:, 0], f["samples"][:, 1] #dset_samples[:, 0], dset_ = f['flags'][...]
    positions = f['positions'][...]
    flat_ids, flat_bps = f['thread_targets'][:, :-1], f['thread_targets'][:, -1]
    # flat_ids, flat_bps = f['thread_targets'][:, 0], f['thread_targets'][:, 1]
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
        Reading in haps and maps file for Li-Stephens
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

def parse_demography(demography):
    d = pd.read_table(demography, delim_whitespace=True, header=None)
    return list(d[0]), list(d[1])
