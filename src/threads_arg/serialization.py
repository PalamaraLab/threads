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

import h5py
import numpy as np
from datetime import datetime

from threads_arg import ThreadingInstructions


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

    if variant_metadata is not None and num_sites:
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
    return ThreadingInstructions(
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
