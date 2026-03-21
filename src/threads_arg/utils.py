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
    phys_list, cm_list, chr_list = [], [], []
    with open(map_file) as f:
        header = f.readline().strip().split()
        pos_idx = header.index('pos')
        chr_idx = header.index('chr')
        cm_idx = header.index('cM')
        for line in f:
            fields = line.strip().split()
            phys_list.append(float(fields[pos_idx]))
            chr_list.append(str(fields[chr_idx]))
            cm_list.append(float(fields[cm_idx]))
    cm_pos = np.array(cm_list, dtype=np.float64)
    phys_pos = np.array(phys_list, dtype=np.float64)
    chromosomes = np.unique(chr_list)

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
    if not pgen_file.endswith("pgen"):
        raise ValueError(f"Cannot find .pvar or .bim files, {pgen_file} does not end with 'pgen'.")
    pvar = pgen_file.rstrip("pgen") + "pvar"
    bim  = pgen_file.rstrip("pgen") + "bim"
    physical_positions = None
    if os.path.isfile(bim):
        pos = []
        with open(bim) as f:
            for line in f:
                if line.startswith('#'): continue
                pos.append(float(line.split()[3]))
        physical_positions = np.array(pos, dtype=np.float64)
    elif os.path.isfile(pvar):
        pos = []
        with open(pvar) as f:
            for line in f:
                if line.startswith('#'): continue
                pos.append(float(line.split()[1]))
        physical_positions = np.array(pos, dtype=np.float64)
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


class VariantMetadata:
    """Lightweight replacement for pandas DataFrame for variant metadata."""
    __slots__ = ('_data', '_len')

    def __init__(self, data):
        self._data = data  # dict of numpy arrays
        self._len = len(next(iter(data.values())))

    def __getitem__(self, key):
        if isinstance(key, str):
            return self._data[key]
        # Boolean mask
        return VariantMetadata({k: v[key] for k, v in self._data.items()})

    def __len__(self):
        return self._len

    @property
    def columns(self):
        return list(self._data.keys())

    @property
    def shape(self):
        return (self._len,)


def read_variant_metadata(pgen):
    """
    Attempt to read variant metadata in vcf style:
    CHR, POS, ID, REF, ALT, QUAL, FILTER
    """
    pvar = pgen.replace("pgen", "pvar")
    bim = pgen.replace("pgen", "bim")
    if os.path.isfile(bim):
        chrom, pos, vid, ref, alt = [], [], [], [], []
        with open(bim) as f:
            for line in f:
                if line.startswith('#'):
                    continue
                fields = line.strip().split()
                chrom.append(fields[0])
                vid.append(fields[1])
                pos.append(fields[3])
                alt.append(fields[4])
                ref.append(fields[5])
        return VariantMetadata({
            "CHROM": np.array(chrom), "POS": np.array(pos),
            "ID": np.array(vid), "REF": np.array(ref), "ALT": np.array(alt),
            "QUAL": np.full(len(pos), "."), "FILTER": np.full(len(pos), "PASS"),
        })
    elif os.path.isfile(pvar):
        # Parse header to find column indices
        header = None
        header_line_count = 0
        with open(pvar) as f:
            for line in f:
                if line.startswith("##"):
                    header_line_count += 1
                    continue
                if line.startswith("#CHROM"):
                    header = line.strip().lstrip('#').split()
                    header_line_count += 1
                    break
        if header is None:
            raise RuntimeError(f"Invalid .pvar file {pvar}")

        col_idx = {name: i for i, name in enumerate(header)}
        chrom, pos, vid, ref, alt, qual, filt = [], [], [], [], [], [], []
        has_filter = "FILTER" in col_idx
        has_qual = "QUAL" in col_idx
        with open(pvar) as f:
            for _ in range(header_line_count):
                next(f)
            for line in f:
                fields = line.strip().split()
                chrom.append(fields[col_idx["CHROM"]])
                pos.append(fields[col_idx["POS"]])
                vid.append(fields[col_idx["ID"]])
                ref.append(fields[col_idx["REF"]])
                alt.append(fields[col_idx["ALT"]])
                filt.append(fields[col_idx["FILTER"]] if has_filter else "PASS")
                qual.append(fields[col_idx["QUAL"]] if has_qual else ".")
        return VariantMetadata({
            "CHROM": np.array(chrom), "POS": np.array(pos),
            "ID": np.array(vid), "REF": np.array(ref), "ALT": np.array(alt),
            "QUAL": np.array(qual), "FILTER": np.array(filt),
        })
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
        with open(psam, "r") as f:
            header_line = f.readline().strip()
            header = header_line.split()
            # Find the IID column
            if "IID" in header:
                iid_idx = header.index("IID")
            elif "#IID" in header:
                iid_idx = header.index("#IID")
            else:
                # No recognized header, treat as fam-like (second column)
                f.seek(0)
                return [l.split()[1] for l in f]
            names = []
            for line in f:
                if line.startswith('#'):
                    continue
                fields = line.strip().split()
                if fields:
                    names.append(fields[iid_idx])
            return names
    else:
        raise RuntimeError(f"Can't find {fam} or {psam}")


def parse_demography(demography):
    times, sizes = [], []
    with open(demography) as f:
        for line in f:
            fields = line.strip().split()
            if len(fields) >= 2:
                times.append(float(fields[0]))
                sizes.append(float(fields[1]))
    return times, sizes


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


def read_all_genotypes(pgen):
    """Read all genotypes from a pgen file into a single (n_sites, n_haps) int32 array."""
    reader = pgenlib.PgenReader(pgen.encode())
    num_samples = reader.get_raw_sample_ct()
    num_sites = reader.get_variant_ct()
    n_haps = 2 * num_samples

    alleles_out = np.empty((num_sites, n_haps), dtype=np.int32)
    phasepresent_out = np.empty((num_sites, num_samples), dtype=np.uint8)

    BATCH_SIZE = max(1, int(4e7 // n_haps))
    for b_start in range(0, num_sites, BATCH_SIZE):
        b_end = min(num_sites, b_start + BATCH_SIZE)
        reader.read_alleles_and_phasepresent_range(
            b_start, b_end,
            alleles_out[b_start:b_end],
            phasepresent_out[b_start:b_end])

    if np.any(phasepresent_out == 0):
        raise RuntimeError("Unphased variants are currently not supported.")

    return alleles_out


def pgen_chunk_iterator(pgen, chunk_size=None):
    """Yield (n_sites_chunk, n_haps) int32 arrays from a pgen file.

    Each chunk is a contiguous block of sites read from disk. The chunks
    are yielded in order and together cover all sites exactly once. The
    buffer is reused across yields for efficiency, so callers must not
    hold references to previous chunks.
    """
    reader = pgenlib.PgenReader(pgen.encode())
    num_samples = reader.get_raw_sample_ct()
    num_sites = reader.get_variant_ct()
    n_haps = 2 * num_samples

    if chunk_size is None:
        chunk_size = max(1, int(4e7 // n_haps))

    alleles_buf = np.empty((chunk_size, n_haps), dtype=np.int32)
    phase_buf = np.empty((chunk_size, num_samples), dtype=np.uint8)

    for b_start in range(0, num_sites, chunk_size):
        b_end = min(num_sites, b_start + chunk_size)
        g_size = b_end - b_start
        out = alleles_buf[:g_size]
        phase_out = phase_buf[:g_size]
        reader.read_alleles_and_phasepresent_range(b_start, b_end, out, phase_out)
        if np.any(phase_out == 0):
            raise RuntimeError("Unphased variants are currently not supported.")
        yield out


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
