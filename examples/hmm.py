from TDPBWT import DPBWT
import numpy as np
import pandas as pd
import msprime
import matplotlib.pyplot as plt

# haploid!
demography = "/well/palamara/users/awo066/TDPBWT/experiments/threading_accuracy/resources/CEU.demo"
SIM_LENGTH = 1e6
RHO = 1.2e-8
MU = 1.65e-8
seed = 130222

df = pd.read_table(demography, header=None)
df.columns  = ['GEN', 'NE']
demography = msprime.Demography()
# NOTE: these initial sizes get overwritten anyways...
demography.add_population(name="A", initial_size=1e4)
for t, ne in zip(df.GEN.values, df.NE.values):
    # Divide by 2 because msprime wants diploids but demography is haploid
    demography.add_population_parameters_change(t, initial_size=ne / 2, population="A")

ts = msprime.sim_ancestry(
    samples={"A": 100},  # diploid!
    demography=demography,
    sequence_length=SIM_LENGTH,
    recombination_rate=RHO,
    random_seed=seed)
mts = msprime.sim_mutations(ts, rate=MU, random_seed=seed)
G = mts.genotype_matrix().transpose()

pos = np.array([s.position for s in mts.sites()])
assert len(pos) == G.shape[1]
cm = RHO * 100 * pos
bwt = DPBWT(pos, cm, MU, list(df.NE), list(df.GEN))

for g in G[:-1]:
    bwt.insert(g)

obs = [1] * 50 + [0] * 500
obs2 = [1] * 10 + [0] * 500 + [1] * 50
print(bwt.hmm.breakpoints(obs, 5))
print(bwt.hmm.breakpoints(obs2, 25))
import pdb
pdb.set_trace()