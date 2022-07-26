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
bw1 = bwt.fastLS([1, 1, 1, 0, 0, 1, 0, 0, 1])
print(bw1)
assert bw1 == [(0, 2)]

print(bwt.thread([1, 1, 1, 0, 0, 1, 0, 0, 1]))

# then something about deletion