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

# print(bwt.longest_prefix([0, 1, 1, 0, 0, 1, 0, 0, 1]))
# print(bwt.longest_prefix([1, 1, 1, 0, 0, 1, 0, 0, 1]))
bwt.print_sorting()
bw1 = bwt.fastLS([1, 1, 1, 0, 0, 1, 0, 0, 1])
print(bw1)
assert bw1 == [(0, 2)]

print(bwt.thread([1, 1, 1, 0, 0, 1, 0, 0, 1]))


# Test 2
print("TEST 2")
pos2 = np.array([1, 2, 3, 4, 5, 6, 7, 8])
bwt = DPBWT(pos2, 100 * rho * pos2, mu, Ne=1)
bwt.insert([0] * 8)
bwt.insert([1] * 8)
bw2a = bwt.fastLS([1, 1, 1, 1, 0, 0, 0, 0])
print(bw2a)
assert bw2a == [(0, 1), (4, 0)]  # [1, 1, 1, 1, 0, 0, 0, 0]

bw2b = bwt.fastLS([1, 1, 0, 0, 0, 0, 1, 1])
print(bw2b)
assert bw2b == [(0, 1), (2, 0), (6, 1)] # [1, 1, 0, 0, 0, 0, 1, 1]

bw2c = bwt.fastLS([1, 0, 1, 0, 1, 0, 1, 1])
print(bw2c)
assert bw2c == [(0, 1)]#[1, 1, 1, 1, 1, 1, 1, 1]

# Test 3
print("\nTEST 3, in which extension by divergence fails")
pos3 = np.array([1, 2, 3, 4, 5, 6, 7, 8])
bwt = DPBWT(pos3, 100 * rho * pos2, mu, Ne=1)
bwt.insert([0, 0, 0, 0, 1, 0, 0, 0])
bwt.insert([0, 0, 0, 0, 1, 0, 0, 0])
bwt.insert([1, 1, 1, 1, 0, 0, 0, 0])
bwt.insert([1, 1, 1, 1, 0, 0, 0, 0])
bw3a = bwt.fastLS([1, 1, 1, 1, 1, 0, 0, 0])
assert bw3a == [(0, 3)]#[3, 3, 3, 3, 3, 3, 3, 3]
print(bw3a)
