from src import DPBWT
bwt = DPBWT([1, 2, 3, 4, 5, 6, 7, 8, 9], [1., 2., 3., 4., 5., 6., 7., 8., 9.], 0.5)
print(bwt)
bwt.insert([0, 1, 0, 0, 0, 1, 0, 0, 1])
bwt.insert([0, 1, 1, 1, 0, 1, 1, 0, 1])
bwt.insert([1, 1, 1, 0, 0, 1, 0, 1, 1])
bwt.insert([1, 0, 0, 1, 0, 1, 1, 0, 1])

print(bwt.longest_prefix([0, 1, 1, 0, 0, 1, 0, 0, 1]))
print(bwt.longest_prefix([1, 1, 1, 0, 0, 1, 0, 0, 1]))
bwt.print_sorting()
bw0 = bwt.fastLS([1, 1, 1, 0, 0, 1, 0, 0, 1])
print(bw0)
assert bw0 == [2] * 9

bwt = DPBWT([1, 2, 3, 4, 5, 6, 7, 8], [1., 2., 3., 4., 5., 6., 7., 8.], 0.5)
bwt.insert([0] * 8)
bwt.insert([1] * 8)
# bwt.print_sorting()
bw1 = bwt.fastLS([1, 1, 1, 1, 0, 0, 0, 0])
print(bw1)
assert bw1 == [1, 1, 1, 1, 0, 0, 0, 0]

bw2 = bwt.fastLS([1, 1, 0, 0, 0, 0, 1, 1])
print(bw2)
assert bw2 == [1, 1, 0, 0, 0, 0, 1, 1]

bw3 = bwt.fastLS([1, 0, 1, 0, 1, 0, 1, 1])
print(bw3)
assert bw3 == [1, 1, 1, 1, 1, 1, 1, 1]
