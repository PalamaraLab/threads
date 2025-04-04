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
from dataclasses import dataclass
from datetime import datetime
import sys

from threads_arg import ThreadingInstructions
from .utils import split_list

@dataclass
class InstructionsData:
    starts: list
    tmrcas: list
    targets: list
    mismatches: list
    positions: list
    region_start: int
    region_end: int


def serialize_instructions(instructions, out, variant_metadata=None, allele_ages=None, sample_names=None):
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

    if variant_metadata is not None:
        # If not none, it is a pandas dataframe with columns
        # CHR, POS, ID, REF, ALT, QUAL, FILTER
        min_pos = min(positions)
        max_pos = max(positions)
        dset_variant_metadata = f.create_dataset("variant_metadata", (num_sites, 7), dtype=h5py.string_dtype(encoding='utf-8'), compression='gzip',
                                    compression_opts=compression_opts)
        variant_metadata = variant_metadata[(variant_metadata["POS"] >= min_pos) & (variant_metadata["POS"] <= max_pos)]
        assert variant_metadata.shape[0] == len(positions)
        assert np.all(np.array(variant_metadata["POS"], dtype=int) == np.array(positions))

        dset_variant_metadata[:, 0] = variant_metadata["CHROM"].astype(str)
        dset_variant_metadata[:, 1] = variant_metadata["POS"].astype(str)
        dset_variant_metadata[:, 2] = variant_metadata["ID"].astype(str)
        dset_variant_metadata[:, 3] = variant_metadata["REF"].astype(str)
        dset_variant_metadata[:, 4] = variant_metadata["ALT"].astype(str)
        dset_variant_metadata[:, 5] = variant_metadata["QUAL"].astype(str)
        dset_variant_metadata[:, 6] = variant_metadata["FILTER"].astype(str)

    if allele_ages is not None:
        assert len(allele_ages) == len(positions)
        dset_allele_ages = f.create_dataset("allele_ages", (num_sites, ), dtype=np.double, compression='gzip',
                                    compression_opts=compression_opts)
        dset_allele_ages[:] = allele_ages

    if sample_names is not None:
        assert len(sample_names) == num_threads // 2
        num_diploids = len(sample_names)
        dset_sample_names = f.create_dataset("sample_names", (num_diploids,), dtype=h5py.string_dtype(encoding='utf-8'), compression='gzip',
                                    compression_opts=compression_opts)
        dset_sample_names[:] = sample_names
    f.close()


def load_instructions(threads):
    """
    Create ThreadingInstructions object from a source .threads file
    """
    inst_data = load_instructions_data(threads)
    return instructions_from_data(inst_data)


def instructions_from_data(instructions_data):
    return ThreadingInstructions(
        instructions_data.starts,
        instructions_data.tmrcas,
        instructions_data.targets,
        instructions_data.mismatches,
        instructions_data.positions,
        instructions_data.region_start,
        instructions_data.region_end
    )


def load_instructions_data_batched(threads, num_batches):
    """
    Load InstructionsData from a .threads file and split into num_batches that
    may be later processed on separate CPUs.

    This returns an array of InstructionsData that may safely passed to new
    processes. The actual use of the C++ API should be done inside each process,
    for example use instructions_from_data on each batch of data and then
    GenotypeIterator and AgeEstimator can be used per-process.
    """
    inst_data = load_instructions_data(threads)

    # Split list will generate cl-open blocks of positions, e.g. [1,2,3,4,5,6]
    # over 2 batches yields batches with positions [1,2,3], [4,5,6], i.e. each
    # batch is complete which is why the in_range() check is inclusive.
    batched_positions = split_list(inst_data.positions, num_batches)

    batched_instructions_data = []
    for bpos in batched_positions:
        range_start = bpos[0]
        range_end = bpos[-1]
        range_starts = []
        range_tmrcas = []
        range_targets = []
        range_mismatches = []

        # Position in range check is inclusive, split_list results do not overlap
        def in_range(pos):
            return (pos >= range_start) and (pos <= range_end)

        # Make positions for this batch as subset of those in range
        range_start_idx = inst_data.positions.index(range_start)
        range_end_idx = inst_data.positions.index(range_end)
        range_positions = inst_data.positions[range_start_idx:range_end_idx + 1]

        # Assumes that tmrcas and targets same size as starts
        for inst_idx in range(len(inst_data.starts)):
            sub_starts = []
            sub_tmrcas = []
            sub_targets = []
            sub_mismatches = []

            # Create sub-vectors based on range
            num_starts = len(inst_data.starts[inst_idx])
            for pos_idx in range(num_starts):
                seg_start = inst_data.starts[inst_idx][pos_idx]
                # FIXME neater way to handle last segment rather than sys.maxsize
                seg_end = sys.maxsize if pos_idx == num_starts - 1 else inst_data.starts[inst_idx][pos_idx + 1]
                if (range_start <= seg_end) and (range_end >= seg_start):
                    sub_starts.append(seg_start)
                    sub_tmrcas.append(inst_data.tmrcas[inst_idx][pos_idx])
                    sub_targets.append(inst_data.targets[inst_idx][pos_idx])

            # Recompute mismatches based on start indexes within range
            for mismatch in inst_data.mismatches[inst_idx]:
                pos = inst_data.positions[mismatch]
                if in_range(pos):
                    sub_mismatches.append(mismatch - range_start_idx)

            range_starts.append(sub_starts)
            range_tmrcas.append(sub_tmrcas)
            range_targets.append(sub_targets)
            range_mismatches.append(sub_mismatches)

        batched_instructions_data.append(
            InstructionsData(
                range_starts,
                range_tmrcas,
                range_targets,
                range_mismatches,
                range_positions,
                range_start,
                range_end
            )
        )

    return batched_instructions_data


def load_instructions_data(threads):
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

    positions = positions.astype(int).tolist()
    return InstructionsData(
        starts,
        tmrcas,
        targets,
        mismatches,
        positions,
        region_start,
        region_end
    )


def load_metadata(threads):
    f = h5py.File(threads, "r")
    import pandas as pd
    return pd.DataFrame(f["variant_metadata"][:], columns=["CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER"])

def load_sample_names(threads):
    f = h5py.File(threads, "r")
    return f["sample_names"][:]
