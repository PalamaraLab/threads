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
import numpy as np

from threads_arg import AgeEstimator, GenotypeIterator
from .serialization import load_instructions
from .utils import timer_block # FIXME remove

def estimate_allele_ages(threads, out):
    logging.info("Starting allele age estimation with the following parameters:")
    logging.info(f"threads:  {threads}")
    logging.info(f"out:      {out}")

    with timer_block("estimate_allele_ages"):
        # Read threading instructions
        instructions = load_instructions(threads)

        # Estimate ages
        gt_it = GenotypeIterator(instructions)
        age_estimator = AgeEstimator(instructions)
        while gt_it.has_next_genotype():
            g = np.array(gt_it.next_genotype())
            age_estimator.process_site(g)
        allele_age_estimates = age_estimator.get_inferred_ages()

        # Temporary snp ids until #45 is resolved
        snp_ids = [f"snp_{i}" for i in range(len(instructions.positions))]

        # Write results to file
        with open(out, "w") as outfile:
            for snp_id, pos, allele_age in zip(snp_ids, instructions.positions, allele_age_estimates):
                outfile.write(f"{snp_id}\t{pos}\t{allele_age}\n")
