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
import numpy as np
from tqdm import tqdm

from threads_arg import AgeEstimator, GenotypeIterator
from .serialization import load_instructions
from .utils import timer_block, default_process_count, split_list

logger = logging.getLogger(__name__)


def _nth_batch_worker(instructions, result_idx, allele_ages_results):
    # Estimate ages on this instruction batch
    gt_it = GenotypeIterator(instructions)
    age_estimator = AgeEstimator(instructions)
    while gt_it.has_next_genotype():
        g = np.array(gt_it.next_genotype())
        age_estimator.process_site(g)
    allele_age_estimates = age_estimator.get_inferred_ages()

    # Index result so full data can be reconstructed in order
    allele_ages_results[result_idx] = allele_age_estimates


def _nth_batch_worker_star(args):
    return _nth_batch_worker(*args)

def estimate_ages(instructions, num_batches, num_threads):
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
        # Process-safe dict so batch results can be reconstructed in order
        manager = multiprocessing.Manager()
        allele_ages_results = manager.dict()

        # Create arguments for each job in process pool
        jobs_args = [(range_inst, i, allele_ages_results)
                        for i, range_inst in enumerate(batched_instructions)]

        with multiprocessing.Pool(processes=num_processors) as pool:
            # To use tqdm with a pool, use imap with shim method to unpack args.
            # Note the enclosing unassigned list() call is necessary. Otherwise
            # no value is retrieved and the process is ignored/dropped.
            list(tqdm(
                pool.imap(_nth_batch_worker_star, jobs_args),
                total=len(jobs_args)
            ))

    # Collect batched estimates into single list in index sort order
    allele_age_estimates = []
    for i in range(len(allele_ages_results)):
        allele_age_estimates += allele_ages_results[i]
    return allele_age_estimates

def estimate_allele_ages(threads, out, num_threads):
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
