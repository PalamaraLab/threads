"""Complete example of sim through tdpbwt arg-inference complete with cycling"""
from TDPBWT import DPBWT
import numpy as np
import pandas as pd
import msprime
import arg_needle_lib

### SIMULATION
demography = "/well/palamara/users/awo066/TDPBWT/experiments/threading_accuracy/resources/CEU.demo" # haploid!
SIM_LENGTH = 1e4
RHO = 1.2e-8
MU = 1.65e-8
seed = 130222
N = 11 # diploid!
N_cycle = 100 # haploid!

df = pd.read_table(demography, header=None)
df.columns  = ['GEN', 'NE']
demography = msprime.Demography()
# NOTE: these initial sizes get overwritten anyways...
demography.add_population(name="A", initial_size=1e4)
for t, ne in zip(df.GEN.values, df.NE.values):
    # Divide by 2 because msprime wants diploids but demography is haploid
    demography.add_population_parameters_change(t, initial_size=ne / 2, population="A")

ts = msprime.sim_ancestry(
    samples={"A": N},  # diploid!
    demography=demography,
    sequence_length=SIM_LENGTH,
    recombination_rate=RHO,
    random_seed=seed)
mts = msprime.sim_mutations(ts, rate=MU, random_seed=seed)
G = mts.genotype_matrix().transpose()
assert G.shape[0] == 2 * N
pos = np.array([s.position for s in mts.sites()])

### ARG INFERENCE
cm = RHO * 100 * pos
bwt = DPBWT(pos, cm, MU, list(df.NE), list(df.GEN))

threading_instructions = []
arg = arg_needle_lib.ARG(0, pos[-1] - pos[0], reserved_samples=2 * N)
arg.set_offset(int(pos[0]))

# How this should work:
for i, g in enumerate(G):
    arg.add_sample(i, str(i))
    if i == 0:
        bwt.insert(g)
    else:
        section_starts, thread_ids, thread_heights = bwt.thread(g)
        arg.thread_sample([s - arg.offset for s in section_starts], thread_ids, thread_heights)

print('wowza, done building first-pass arg')
# Cycling pass
for i, g in enumerate(G[:N_cycle]):
    bwt.delete_ID(i)
    arg.lazy_delete(i)
    arg.add_sample(i, str(i))
    section_starts, thread_ids, thread_heights = bwt.thread(i, g)
    arg.thread_sample([s - arg.offset for s in section_starts], thread_ids, thread_heights)
print('wowza, done cycling')
