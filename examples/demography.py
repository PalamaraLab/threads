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
coal_distr = [] 
coal_x = np.arange(0.05, 2.51, 0.01)
for t in coal_x:
    coal_distr.append(bwt.demography.std_to_gen(t))

plt.plot(coal_x, coal_distr)
plt.yscale('log')
plt.savefig('example_data/coal.png')
for g in G[:-1]:
    bwt.insert(g)

thread = bwt.thread(G[-1])
print(thread)

plt.cla()
Ns = np.arange(2, 1_000)
sims = []
for n in Ns:
    ts = msprime.sim_ancestry(
    samples={"A": n},
    ploidy=1,
    demography=demography)
    tree = ts.first()
    sims.append(tree.time(tree.parent(0)))
exp_length = [bwt.demography.expected_branch_length(i) for i in Ns]
exp_length_2 = [bwt.demography.std_to_gen(2 / i) for i in Ns]
exp_length_3 = [bwt.demography.std_to_gen(1 / i) for i in Ns]
plt.plot(Ns, exp_length, label='exp')
plt.scatter(Ns, sims, label='sim')
plt.plot(Ns, exp_length_2, label='stdgen')
plt.plot(Ns, exp_length_3, label='stdgen1')
plt.yscale('log')
plt.legend()
plt.savefig("example_data/exp_lengths.png")
