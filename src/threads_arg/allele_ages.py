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

import logging
import multiprocessing
import numpy as np

from threads_arg import AgeEstimator, GenotypeIterator
from .serialization import load_instructions_data_batched, instructions_from_data
from .utils import timer_block, default_process_count

logger = logging.getLogger(__name__)


def _allele_ages_worker(instructions_data, result_idx, allele_ages_results):
    with timer_block(f"processing instructions for batch {result_idx}", print_start=False):
        instructions = instructions_from_data(instructions_data)

        # Estimate ages on this instruction batch
        gt_it = GenotypeIterator(instructions)
        age_estimator = AgeEstimator(instructions)
        while gt_it.has_next_genotype():
            g = np.array(gt_it.next_genotype())
            age_estimator.process_site(g)
        allele_age_estimates = age_estimator.get_inferred_ages()

        # Index result so full data can be reconstructed in order
        allele_ages_results[result_idx] = allele_age_estimates


def estimate_allele_ages(threads, out, num_batches):
    # If not specified use 75% of CPUs available; performance diminishes at 100%
    if not num_batches:
        num_batches = int(default_process_count() * 0.75)

    logging.info("Starting allele age estimation with the following parameters:")
    logging.info(f"threads:     {threads}")
    logging.info(f"out:         {out}")
    logging.info(f"num_batches: {num_batches}")

    with timer_block(f"loading {threads} into {num_batches} batches", print_start=False):
        batched_instructions_data = load_instructions_data_batched(threads, num_batches)
        positions = []
        for batch in batched_instructions_data:
            positions += batch.positions

    with timer_block(f"estimating allele ages ({num_batches} CPUs)"):
        # Process-safe dict so batch results can be reconstructed in order
        manager = multiprocessing.Manager()
        allele_ages_results = manager.dict()
        jobs = []
        for i in range(num_batches):
            p = multiprocessing.Process(
                target=_allele_ages_worker,
                args=(batched_instructions_data[i], i, allele_ages_results)
            )
            jobs.append(p)
            p.start()

        for p in jobs:
            p.join()

    logger.info(f"Writing results to {out}...")

    # Collect batched estimates into single list in index sort order
    allele_age_estimates = []
    for i in range(len(allele_ages_results)):
        allele_age_estimates += allele_ages_results[i]

    # Temporary snp ids until #45 is resolved
    snp_ids = [f"snp_{i}" for i in range(len(positions))]

    # Write results to file
    with open(out, "w") as outfile:
        for snp_id, pos, allele_age in zip(snp_ids, positions, allele_age_estimates):
            outfile.write(f"{snp_id}\t{pos}\t{allele_age}\n")
