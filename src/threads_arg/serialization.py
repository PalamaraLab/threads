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

import numpy as np

from threads_arg import (
    serialize_threads as _serialize_threads,
    deserialize_threads as _deserialize_threads,
    read_threads_metadata as _read_threads_metadata,
    read_threads_sample_names as _read_threads_sample_names,
)


def serialize_instructions(instructions, out, variant_metadata=None, allele_ages=None, sample_names=None):
    """Serialize ThreadingInstructions to a .threads file.

    Args:
        instructions: ThreadingInstructions object.
        out: Output file path (without extension).
        variant_metadata: Optional VariantMetadata with CHROM/POS/ID/REF/ALT/QUAL/FILTER.
        allele_ages: Optional list of allele age estimates per site.
        sample_names: Optional list of sample names (one per diploid individual).
    """
    positions = instructions.positions

    metadata_cols = []
    if variant_metadata is not None and instructions.num_sites:
        min_pos = min(positions)
        max_pos = max(positions)
        pos_int = np.array(variant_metadata["POS"], dtype=int)
        mask = (pos_int >= min_pos) & (pos_int <= max_pos)
        variant_metadata = variant_metadata[mask]
        assert variant_metadata.shape[0] == len(positions)
        assert np.all(np.array(variant_metadata["POS"], dtype=int) == np.array(positions))
        for col in ["CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER"]:
            metadata_cols.append([str(x) for x in variant_metadata[col]])

    ages = list(allele_ages) if allele_ages is not None else []
    names = [str(x) for x in sample_names] if sample_names is not None else []

    _serialize_threads(out, instructions, metadata_cols, ages, names)


def load_instructions(threads):
    """Load ThreadingInstructions from a .threads file."""
    return _deserialize_threads(threads)


def load_metadata(threads):
    """Load variant metadata (CHROM, POS, ID, REF, ALT, QUAL, FILTER) from a .threads file."""
    from .utils import VariantMetadata
    cols = _read_threads_metadata(threads)
    columns = ["CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER"]
    return VariantMetadata({col: np.array(data) for col, data in zip(columns, cols)})


def load_sample_names(threads):
    """Load sample names from a .threads file."""
    return _read_threads_sample_names(threads)
