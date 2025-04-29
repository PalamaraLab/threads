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

import os
import numpy as np
import pandas as pd
import logging
import pgenlib
import re
import time
import warnings

from contextlib import contextmanager
from typing import Tuple, Union

logger = logging.getLogger(__name__)

def read_map_file(map_file, expected_chromosome=None) -> Tuple[np.ndarray, np.ndarray, str]:
    """
    Reading in map file for Li-Stephens using genetic maps in the SHAPEIT format
    """
    maps = pd.read_table(map_file, sep=r"\s+")
    cm_pos = maps.cM.values.astype(np.float64)
    phys_pos = maps.pos.values.astype(np.float64)
    chromosomes = np.unique(maps.chr.values.astype(str))

    # Currently we only allow for processing one chromosome at a time
    if len(chromosomes) > 1:
        raise RuntimeError(f"Found multiple chromosomes in {map_file}. Threads limits to one.")

    if expected_chromosome:
        if str(expected_chromosome) != chromosomes[0]:
            raise RuntimeError(f"Expected chromosome {expected_chromosome} not found in {map_file}")

    # Add an epsilon value to avoid div zero errors. A few decimal places is
    # enough to distinguish each cM, so 1e-5 does not skew results.
    for i in range(1, len(cm_pos)):
        if cm_pos[i] <= cm_pos[i-1]:
            cm_pos[i] = cm_pos[i-1] + 1e-5
    return phys_pos, cm_pos, chromosomes[0]


def _read_pgen_physical_positions(pgen_file):
    pvar = pgen_file.replace("pgen", "pvar")
    bim = pgen_file.replace("pgen", "bim")
    physical_positions = None
    if os.path.isfile(bim):
        physical_positions = np.array(pd.read_table(bim, sep="\\s+", header=None, comment='#')[3]).astype(np.float64)
    elif os.path.isfile(pvar):
        physical_positions = np.array(pd.read_table(pvar, sep="\\s+", header=None, comment='#')[1]).astype(np.float64)
    else:
        raise RuntimeError(f"Can't find {bim} or {pvar} for {pgen_file}")

    return physical_positions


def make_recombination_from_map_and_pgen(map_file, pgen_file, expected_chrom):
    """
    Interpolate pgen to cM positons from map file in SHAPEIT and pgen variants
    """
    phys_pos, cm_pos, _chrom = read_map_file(map_file, expected_chrom)

    physical_positions = _read_pgen_physical_positions(pgen_file)
    cm_out = np.interp(physical_positions, phys_pos, cm_pos)

    if physical_positions.max() > phys_pos.max() or physical_positions.min() < phys_pos.min():
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


def read_variant_metadata(pgen):
    """
    Attempt to read variant metadata in vcf style:
    CHR, POS, ID, REF, ALT, QUAL, FILTER
    """
    pvar = pgen.replace("pgen", "pvar")
    bim = pgen.replace("pgen", "bim")
    if os.path.isfile(bim):
        bim_df = pd.read_table(bim, names=["CHROM", "ID", "CM", "POS", "ALT", "REF"])
        out_df = bim_df[["CHROM", "POS", "ID", "REF", "ALT"]]
        out_df["FILTER"] = bim_df["FILTER"] if "FILTER" in out_df.columns else "PASS"
        out_df["QUAL"] = bim_df["QUAL"] if "QUAL" in out_df.columns else "."
        return out_df
    elif os.path.isfile(pvar):
        header = None
        with open(pvar, "r") as pvarfile:
            for line in pvarfile:
                if line.startswith("##"):
                    continue
                if line.startswith("#CHROM"):
                    header = line.strip().split()
                    break
            if header is None:
                raise RuntimeError(f"Invalid .pvar file {pvar}")
        pvar_df = pd.read_table(pvar, comment="#", header=None, names=header, sep=r"\s+").rename({"#CHROM": "CHROM"}, axis=1)

        pd.options.mode.chained_assignment = None
        out_df = pvar_df[["CHROM", "POS", "ID", "REF", "ALT"]]
        out_df["FILTER"] = pvar_df["FILTER"].copy() if "FILTER" in out_df.columns else "PASS"
        out_df["QUAL"] = pvar_df["QUAL"].copy() if "QUAL" in out_df.columns else "."
        return out_df
    else:
        raise RuntimeError(f"Can't find {bim} or {pvar}")

def make_constant_recombination_from_pgen(pgen_file, rho):
    """
    Read pgen variant file and generate a constant recombination using rho
    """
    physical_positions = _read_pgen_physical_positions(pgen_file)
    cm_out = rho * 100 * physical_positions

    for i in range(1, len(cm_out)):
        if cm_out[i] <= cm_out[i-1]:
            cm_out[i] = cm_out[i-1] + 1e-5
    return cm_out, physical_positions

def read_sample_names(pgen):
    """
    Read the sample names corresponding to the input pgen
    """
    fam = pgen.replace("pgen", "fam")
    psam = pgen.replace("pgen", "psam")
    if os.path.isfile(fam):
        with open(fam, "r") as famfile:
            return [l.split()[1] for l in famfile]

    elif os.path.isfile(psam):
        sam_df = pd.read_table(psam, sep=r"\s+")
        if "IID" in sam_df.columns:
            return sam_df["IID"].astype(str).tolist()
        elif "#IID" in sam_df.columns:
            return sam_df["#IID"].astype(str).tolist()
        else:
            # If no header, default to famfile
            with open(psam, "r") as famfile:
                return [l.split()[1] for l in famfile]
    else:
        raise RuntimeError(f"Can't find {fam} or {psam}")


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


def parse_region_string(region: str) -> Tuple[Union[str, None], int, int]:
    """
    Convert "chr:start-end" or "start-end" string into a (str | None, int, int)
    3-tuple. The format of chr is either "chr[1-22]" or just "[1-22]". If the
    chr is omitted then the returned 3-tuple's first value is None.
    """
    match = re.match(r"((chr)?(\d+):)?(\d+)-(\d+)", region)
    if not match:
        raise RuntimeError(f"Invalid region string '{region}'")

    chr = match[3]
    if chr:
        chr = int(match[3])
        if chr < 1 or chr > 22:
            raise RuntimeError(f"Invalid chromosome {chr} not between 1 and 22")

    start = int(match[4])
    end = int(match[5])
    if end < start:
        raise RuntimeError(f"Invalid range: end {end} less than start {start}")

    return chr, start, end


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


def default_process_count():
    """
    Get the number of CPUs available for multi-processing work

    This tries os.sched_getaffinity first for common case of running on an HPC
    node, e.g. for `srun` with `--ntasks-per-node 4`, this returns 4.
    Some platforms like macOS do not provide this method. In which case, assume
    that this is not an HPC node and just use the host's CPU count directly.
    """
    try:
        return len(os.sched_getaffinity(0))
    except AttributeError:
        return os.cpu_count()
