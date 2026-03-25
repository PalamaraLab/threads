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

import logging
import multiprocessing
from tqdm import tqdm

from threads_arg import estimate_ages as _estimate_ages_cpp
from .serialization import load_instructions
from .utils import timer_block, default_process_count, split_list

logger = logging.getLogger(__name__)


def _batch_worker(instructions):
    return _estimate_ages_cpp(instructions)

def estimate_ages(instructions, num_batches, num_threads):
    """Estimate allele ages from threading instructions using parallel batches.

    Args:
        instructions: ThreadingInstructions object.
        num_batches: Number of genomic sub-ranges to split into (0 = auto).
        num_threads: Maximum number of worker processes.

    Returns:
        List of allele age estimates, one per site.
    """
    # Make sure we don't use more CPUs than requested
    if not num_threads:
        num_threads = 1
    num_processors = min(num_threads, max(1, default_process_count() - 1))

    # Between 2 and 3 batches per processor in the pool seems to be a good
    # default for balancing moving data around and sharing workload.
    if not num_batches:
        num_batches = num_processors * 3

    with timer_block(f"Splitting instructions into {num_batches} batches", print_start=False):
        # Load instructions, and get sub-ranges based on number of batches
        batched_instructions = []
        batch_positions = split_list(instructions.positions, num_batches)
        for bpos in batch_positions:
            range_start = bpos[0]
            range_end = bpos[-1]
            range_instructions = instructions.sub_range(range_start, range_end)
            batched_instructions.append(range_instructions)

    with timer_block(f"Estimating allele ages ({num_processors} CPUs)"):
        with multiprocessing.Pool(processes=num_processors) as pool:
            batch_results = list(tqdm(
                pool.imap(_batch_worker, batched_instructions),
                total=len(batched_instructions)
            ))

    allele_age_estimates = []
    for batch in batch_results:
        allele_age_estimates += batch
    return allele_age_estimates

def estimate_allele_ages(threads, out, num_threads):
    """CLI entry point: estimate allele ages and write to a TSV file."""
    logging.info("Starting allele age estimation with the following parameters:")
    logging.info(f"threads:     {threads}")
    logging.info(f"out:         {out}")
    logging.info(f"num_threads: {num_threads}")
    num_batches = 3 * num_threads
    
    instructions = load_instructions(threads)
    allele_age_estimates = estimate_ages(instructions, num_batches, num_threads)

    # Temporary snp ids until #45 is resolved
    snp_ids = [f"snp_{i}" for i in range(len(instructions.positions))]

    # Write results to file
    logger.info(f"Writing results to {out}...")
    with open(out, "w") as outfile:
        for snp_id, pos, allele_age in zip(snp_ids, instructions.positions, allele_age_estimates):
            outfile.write(f"{snp_id}\t{pos}\t{allele_age}\n")
