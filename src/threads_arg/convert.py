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

import sys
import time
import tszip
import logging
import arg_needle_lib

import numpy as np

from .utils import decompress_threads

logger = logging.getLogger(__name__)


def threads_to_arg(thread_dict, noise=0.0, max_n=None, verify=False, random_seed=0):
    """
    Assemble threading instructions into an ARG
    """
    threading_instructions = thread_dict['threads']
    pos = thread_dict['positions']
    N = max_n if max_n is not None else np.max(thread_dict['samples'].astype(int)) + 1
    logger.info(f"Will thread {N} haploids")
    arg_start, arg_end = thread_dict['arg_range']
    # "+ 2" so we can include mutations on the last site, ARG is end-inclusive, we're not.
    arg = arg_needle_lib.ARG(0, arg_end - arg_start + 2, reserved_samples=N)
    arg.set_offset(int(arg_start))

    rng = np.random.default_rng(seed=random_seed)

    # How this should work:
    for i, t in enumerate(threading_instructions[:N]):
        if i % 1000 == 0:
            logger.info(f"Sequence {i + 1}...")
        arg.add_sample(str(i))
        if i > 0:
            section_starts, thread_ids, thread_heights, _ = t
            if len(thread_ids.shape) == 2:
                thread_ids = thread_ids[:, 0]
            thread_heights += thread_heights * rng.normal(0.0, noise, len(thread_heights))

            # arg will throw exception if there is a collision in heights. In this instance,
            # the caller will increase the amount of noise to offset further and try again.
            arg_starts = [s - arg.offset for s in section_starts]
            try:
                if arg_starts[-1] >= arg.end:
                    arg.thread_sample([s - arg.offset for s in section_starts[:-1]], thread_ids[:-1], thread_heights[:-1])
                else:
                    arg.thread_sample([s - arg.offset for s in section_starts], thread_ids, thread_heights)
            except ValueError:
                import pdb
                pdb.set_trace()
    logger.info(f"Done threading")

    if verify:
        logger.info("Verifying ARG...")
        arg.check_basic()
    return arg


# Implementation is separated from Click entrypoint for use in tests
def threads_convert(threads, argn, tsz, max_n, random_seed, verify):
    """
    Convert input .threads file into .threads or .argn file
    """
    start_time = time.time()
    logger.info(f"Starting Threads-convert with the following parameters:")
    logger.info(f"  threads:     {threads}")
    logger.info(f"  argn:        {argn}")
    logger.info(f"  tsz:         {tsz}")
    logger.info(f"  max_n:       {max_n}")
    logger.info(f"  random_seed: {random_seed}")

    if argn is None and tsz is None:
        logger.info("Nothing to do, quitting.")
        sys.exit(0)
    decompressed_threads = decompress_threads(threads)
    try:
        logger.info("Attempting to convert to arg format...")
        arg = threads_to_arg(decompressed_threads, noise=0.0, max_n=max_n, verify=verify, random_seed=random_seed)
    except:
        # arg_needle_lib does not allow polytomies
        logger.info(f"Conflicting branches (this is expected), retrying with noise=1e-6...")
        try:
            arg = threads_to_arg(decompressed_threads, noise=1e-6, max_n=max_n, verify=verify, random_seed=random_seed)
        except:# tskit.LibraryError:
            logger.info(f"Conflicting branches, retrying with noise=1e-3...")
            arg = threads_to_arg(decompressed_threads, noise=1e-3, max_n=max_n, verify=verify, random_seed=random_seed)
    if argn is not None:
        logger.info(f"Writing to {argn}")
        arg_needle_lib.serialize_arg(arg, argn)
    if tsz is not None:
        logger.info(f"Converting to tree sequence and writing to {tsz}")
        tszip.compress(arg_needle_lib.arg_to_tskit(arg), tsz)
    logger.info(f"Done, in {time.time() - start_time} seconds")
