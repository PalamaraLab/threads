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

from threads_arg import ConsistencyWrapper, AgeEstimator
from .serialization import load_instructions, serialize_instructions
from .utils import iterate_pgen, read_positions_and_ids


def fit_to_data(threads, pgen, region, allele_ages, out):
    logging.info("Starting data consistency routine with the following parameters:")
    logging.info(f"threads:      {threads}")
    logging.info(f"pgen:         {pgen}")
    logging.info(f"region:       {region}")
    logging.info(f"allele_ages:  {allele_ages}")
    logging.info(f"out:          {out}")
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

    region_start = max([arg_start, region_start, pgen_start])
    region_end = min([arg_end, region_end, pgen_end])

    start_idx = np.searchsorted(positions, region_start)
    end_idx = np.searchsorted(positions, region_end, side="right")
    logging.info(f"Will make threading instructions consistent with {end_idx - start_idx} variants")

    logging.info("Starting data-consistency post-processing")
    start_idx = np.searchsorted(positions, region_start)
    end_idx = np.searchsorted(positions, region_end, side="right")

    if allele_ages is None:
        logging.info(f"Inferring allele ages from data")
        age_estimator = AgeEstimator(instructions)
        iterate_pgen(pgen, lambda i, g: age_estimator.process_site(g), start_idx=start_idx, end_idx=end_idx)
        allele_age_estimates = age_estimator.get_inferred_ages()
        assert len(allele_age_estimates) == len(instructions.positions)
    else:
        allele_age_estimates = []
        with open(allele_ages, "r") as agefile:
            for line in agefile:
                agefile_id, agefile_estimate = line.strip().split()[0], line.strip().split()[1]
                if agefile_id not in ids:
                    continue
                allele_age_estimates.append(agefile_estimate)
        try:
            assert len(allele_age_estimates) == end_idx - start_idx
        except AssertionError:
            raise RuntimeError(f"Allele age estimates do not match markers in the region requested.")

    # Start the consistifying
    cw = ConsistencyWrapper(instructions, allele_age_estimates)
    iterate_pgen(pgen, lambda i, g: cw.process_site(g), start_idx=start_idx, end_idx=end_idx)
    consistent_instructions = cw.get_consistent_instructions()
    serialize_instructions(consistent_instructions, out)
