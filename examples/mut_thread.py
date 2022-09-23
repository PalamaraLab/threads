from TDPBWT import DPBWT
import numpy as np
mu = 0.001
rho = 0.001

# Test 1
pos1 = np.array([1, 2, 3, 4, 5, 6, 7, 8, 9])
bwt = DPBWT(pos1, 100 * rho * pos1, mu, Ne=1)
print(bwt)
bwt.insert([0, 1, 0, 0, 0, 1, 0, 0, 1])
bwt.insert([0, 1, 1, 1, 0, 1, 1, 0, 1])
bwt.insert([1, 1, 1, 0, 0, 1, 0, 1, 1])
bwt.insert([1, 0, 0, 1, 0, 1, 1, 0, 1])

bwt.print_sorting()
bw1 = bwt.thread_with_mutations([1, 1, 1, 0, 0, 1, 0, 0, 1])
print(bw1)
# import pdb
# pdb.set_trace()
assert bw1 == ([1], [2], [2.8957528957528953], [8], [False])

# Test 2
print("TEST 2")
pos2 = np.array([1, 2, 3, 4, 5, 6, 7, 8])
bwt = DPBWT(pos2, 100 * rho * pos2, mu, Ne=1)
bwt.insert([0] * 8)
bwt.insert([1] * 8)
_, _, _, het_sites, het_types = bwt.thread_with_mutations([1, 1, 0, 0, 0, 0, 1, 1])
assert len(het_sites) == 0
assert len(het_types) == 0

# Test 3
print("\nTEST 3,")
pos3 = np.array([1, 2, 3, 4, 5, 6, 7, 8])
bwt = DPBWT(pos3, 100 * rho * pos2, mu, Ne=1)
bwt.insert([0, 0, 0, 0, 1, 0, 0, 0])
bwt.insert([0, 0, 0, 0, 1, 0, 0, 0])
bwt.insert([1, 1, 1, 1, 0, 0, 0, 0])
bwt.insert([1, 1, 1, 1, 0, 0, 0, 0])
bw3a = bwt.thread_with_mutations([1, 1, 1, 1, 1, 0, 0, 0])

print(bw3a)
