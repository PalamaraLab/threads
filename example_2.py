from TDPBWT import DPBWT
import numpy as np
import msprime

Ne = 10_000 # (haploid)
SIM_LENGTH = 5e6
RHO = 1.2e-8
MU = 1.65e-8
seed = 130222

ts = msprime.sim_ancestry(
    samples=100,
    population_size=Ne / 2,
    sequence_length=SIM_LENGTH,
    recombination_rate=RHO,
    random_seed=seed)
mts = msprime.sim_mutations(ts, rate=MU, random_seed=seed)

G = mts.genotype_matrix().transpose()
pos = np.array([s.position for s in mts.sites()])
assert len(pos) == G.shape[1]
cm = RHO * 100 * pos

bwt = DPBWT(pos, cm, MU, Ne)
for g in G[:-1]:
    bwt.insert(g)

thread = bwt.thread(G[-1])
import pdb
pdb.set_trace()
