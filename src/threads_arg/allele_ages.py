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

from threads_arg import AgeEstimator
from .serialization import load_instructions
from .utils import iterate_pgen, read_positions_and_ids

def estimate_allele_ages(threads, pgen, region, out):
    logging.info("Starting allele age estimation with the following parameters:")
    logging.info(f"threads:  {threads}")
    logging.info(f"pgen:     {pgen}")
    logging.info(f"region:   {region}")
    logging.info(f"out:      {out}")
    # Read threading instructions
    instructions = load_instructions(threads)

    # Find the intersection of ARG/pgen/requested regions
    if region is not None:
        region_start = int(region.split(":")[-1].split("-")[0])
        region_end = int(region.split(":")[-1].split("-")[1])
    else:
        region_start = -np.inf
        region_end = np.inf
    arg_start = instructions.start
    arg_end = instructions.end

    positions, ids = read_positions_and_ids(pgen)

    pgen_start = positions[0]
    pgen_end = positions[-1] + 1

    allele_age_start = max([arg_start, region_start, pgen_start])
    allele_age_end = min([arg_end, region_end, pgen_end])

    start_idx = np.searchsorted(positions, allele_age_start)
    end_idx = np.searchsorted(positions, allele_age_end, side="right")
    logging.info(f"Will estimate the age of {end_idx - start_idx} variants")

    # Initialize the age estimator
    age_estimator = AgeEstimator(instructions)

    # Do the age estimation
    iterate_pgen(pgen, lambda i, g: age_estimator.process_site(g), start_idx=start_idx, end_idx=end_idx)
    allele_age_estimates = age_estimator.get_inferred_ages()

    # Write results to file
    with open(out, "w") as outfile:
        for allele_age, snp_id in zip(allele_age_estimates, ids[start_idx:end_idx]):
            outfile.write(f"{snp_id}\t{allele_age}\n")

