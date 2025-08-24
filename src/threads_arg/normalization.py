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
import pandas as pd
import math
import numpy as np
import pandas
from threads_arg import ThreadingInstructions
import msprime

logging.basicConfig(
    format='%(asctime)s %(levelname)-8s %(message)s',
    level=logging.INFO,
    datefmt='%Y-%m-%d %H:%M:%S')
logging.getLogger("msprime").setLevel(logging.WARNING)

class Normalizer:
    """
    Minimal simulator based on ARG-Needle Simulator class
    and including the ARG-Needle normalisation routine.
    """
    def __init__(self, demography_file, num_samples):
        self.num_samples = num_samples
        self.demography = self.read_demography(demography_file)

    def read_demography(self, demography_file):
        df = pd.read_table(demography_file, header=None)
        df.columns  = ['GEN', 'NE']
        demography = msprime.Demography()
        # NOTE: these initial sizes get overwritten anyways...
        demography.add_population(name="A", initial_size=1e4)
        for t, ne in zip(df.GEN.values, df.NE.values):
            # Divide by 2 because msprime wants diploids but demography is haploid
            demography.add_population_parameters_change(t, initial_size=ne / 2, population="A")
        return demography

    def simulation(self, length, random_seed=10):
        ts = msprime.sim_ancestry(
            samples={"A": self.num_samples},
            demography=self.demography,
            sequence_length=length,
            recombination_rate=0,
            random_seed=random_seed)
        return ts

    def normalize(self, threading_instructions, num_seeds=1000, start_seed=1):
        num_samples = threading_instructions.num_samples
        assert num_samples == self.num_samples
        heights_to_span = {}
        # could rewrite this using numpy
        positions = threading_instructions.positions

        for start_vec, tmrca_vec in zip(threading_instructions.all_starts(), threading_instructions.all_tmrcas()):
            assert len(start_vec) == len(tmrca_vec)
            for i, (start, tmrca) in enumerate(zip(start_vec, tmrca_vec)):

                try:
                    heights_to_span[tmrca] = start_vec[i + 1] - start
                except IndexError:
                    heights_to_span[tmrca] = positions[-1] + 1 - start

        numpy_stuff = np.array(sorted(heights_to_span.items()))

        cumsum = np.cumsum(numpy_stuff[:, 1])
        quantiles = (cumsum - 0.5 * numpy_stuff[:, 1]) / cumsum[-1]
        del heights_to_span
        del cumsum

        node_times = []
        highest_times = []

        for seed_offset in range(num_seeds):
            seed = seed_offset + start_seed

            # 1e6 doesn't matter here as rho = 0
            simulation = self.simulation(1e3, random_seed=seed)
            highest = 0
            for node in simulation.nodes():
                if node.time > 0:
                    node_times.append(node.time)
                    highest = max(node.time, highest)
            highest_times.append(highest)

        highest_times.sort()
        node_times.sort()

        node_times = np.array(node_times)
        node_times_pad = np.zeros(len(node_times) + 2)
        node_times_pad[1:-1] = np.array(node_times)
        node_times_pad[0] = node_times[0] * 0  # 0 is tunable, as long as it's < 1
        node_times_pad[-1] = node_times[-1] * 1.05  # 1.05 is tunable, as long as it's > 1
        sim_quantiles = np.linspace(0, 1, len(node_times_pad))
        del node_times

        corrected = np.interp(quantiles, sim_quantiles, node_times_pad)
        correction_dict = dict(zip(numpy_stuff[:, 0], corrected))
        del quantiles
        del sim_quantiles
        del node_times_pad
        del numpy_stuff
        del corrected

        new_tmrcas = []
        for tmrca_vec in threading_instructions.all_tmrcas():
            new_tmrca_vec = []
            for tmrca in tmrca_vec:
                new_tmrca_vec.append(correction_dict[tmrca])
            new_tmrcas.append(new_tmrca_vec)

        return ThreadingInstructions(
            threading_instructions.all_starts(),
            new_tmrcas,
            threading_instructions.all_targets(),
            threading_instructions.all_mismatches(),
            positions,
            threading_instructions.start,
            threading_instructions.end
        )
