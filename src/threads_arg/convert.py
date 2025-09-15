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

import sys
import time
import tszip
import logging
import threads_arg
import arg_needle_lib

import numpy as np

from .serialization import load_instructions

logger = logging.getLogger(__name__)


def threads_to_arg(instructions, add_mutations=False, noise=0.0, polytomy_threshold=None, defragment=False):
    """
    Assemble threading instructions into an ARG
    """
    N = instructions.num_samples
    logger.info(f"Will thread {N} haploids")
    arg_start, arg_end = instructions.start, instructions.end
    # "+ 2" so we can include mutations on the last site, ARG is end-inclusive, we're not.
    arg = arg_needle_lib.ARG(0, arg_end - arg_start + 2, reserved_samples=N)
    arg.set_offset(int(arg_start))

    rng = np.random.default_rng(seed=1234)

    # If polytomies are not allowed, fall back to original code that throws an error if matching heights found
    allow_polytomy = polytomy_threshold is not None

    # How this should work:
    for i, (section_starts, thread_ids, thread_heights) in enumerate(zip(instructions.all_starts(), instructions.all_targets(), instructions.all_tmrcas())):
        if i == N:
            break
        arg.add_sample(str(i))
        if i > 0:
            # If polytomy is not allowed, the ARG will throw exception if there is a collision in heights.
            # In this instance, the caller will increase the amount of noise to offset further and try again.
            thread_heights = np.array(thread_heights)
            if not allow_polytomy:
                thread_heights += thread_heights * rng.normal(0.0, noise, len(thread_heights))

            arg_starts = [s - arg.offset for s in section_starts]
            if arg_starts[-1] >= arg.end:
                arg.thread_sample([s - arg.offset for s in section_starts[:-1]], thread_ids[:-1], thread_heights[:-1], polytomy_threshold)
            else:
                arg.thread_sample([s - arg.offset for s in section_starts], thread_ids, thread_heights, polytomy_threshold)
    logger.info(f"Done threading")

    # Matching heights during threading cause stacked edges with same height as parent, these are collapsed to
    # group together at the heighest parent. Node count reported as measure of how much this has changed ARG
    if allow_polytomy:
        logger.info(f"Grouping polytomies within height threshold {polytomy_threshold} (current node count {arg.num_nodes()})...")
        arg.collapse_polytomies(polytomy_threshold)
        logger.info(f"Node count after grouping polytomies {arg.num_nodes()}")

    # Threading can cause adjacent edges that have same parent and child nodes, this is the most
    # simple form of fragmentation and can be fixed by merging adjacent spans if found.
    if defragment:
        logger.info(f"Fixing fragmentation (currently {arg.fragmentation_count()} edges)...")
        arg.fix_fragmentation()

    if add_mutations:
        arg.populate_children_and_roots()
        gt_it = threads_arg.GenotypeIterator(instructions)
        positions = instructions.positions
        i = 0
        while gt_it.has_next_genotype():
            g = gt_it.next_genotype()
            arg_needle_lib.map_genotype_to_ARG(arg, g, positions[i] - arg.offset)
            i += 1

    return arg


# Implementation is separated from Click entrypoint for use in tests
def threads_convert(threads, argn, tsz, add_mutations=False, polytomy_threshold=None, defragment=False):
    """
    Convert input .threads file into .threads or .argn file
    """
    start_time = time.time()
    logger.info(f"Starting Threads-convert with the following parameters:")
    logger.info(f"  threads:       {threads}")
    logger.info(f"  argn:          {argn}")
    logger.info(f"  tsz:           {tsz}")
    logger.info(f"  add_mutations: {add_mutations}")

    if argn is None and tsz is None:
        logger.info("Nothing to do, quitting.")
        sys.exit(0)
    instructions = load_instructions(threads)
    if polytomy_threshold is None:
        try:
            logger.info("Attempting to convert to arg format...")
            arg = threads_to_arg(instructions, add_mutations=add_mutations, noise=0.0, defragment=defragment)
        except:
            # arg_needle_lib does not allow polytomies
            logger.info(f"Conflicting branches (this is expected), retrying with noise=1e-5...")
            try:
                arg = threads_to_arg(instructions, add_mutations=add_mutations, noise=1e-5, defragment=defragment)
            except:# tskit.LibraryError:
                logger.info(f"Conflicting branches, retrying with noise=1e-3...")
                arg = threads_to_arg(instructions, add_mutations=add_mutations, noise=1e-3, defragment=defragment)
    else:
        logger.info(f"Convert to arg format, allowing for polytomies...")
        arg = threads_to_arg(
            instructions,
            add_mutations=add_mutations,
            polytomy_threshold=float(polytomy_threshold),
            defragment=defragment
        )

    if argn is not None:
        logger.info(f"Writing to {argn}")
        arg_needle_lib.serialize_arg(arg, argn)
    if tsz is not None:
        logger.info(f"Converting to tree sequence and writing to {tsz}")
        tszip.compress(arg_needle_lib.arg_to_tskit(arg), tsz)
    logger.info(f"Done, in {time.time() - start_time} seconds")
