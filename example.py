from src import DPBWT
bwt = DPBWT([1, 2, 3, 4, 5, 6, 7, 8, 9], [1., 2., 3., 4., 5., 6., 7., 8., 9.], 0.5)
print(bwt)
bwt.insert([0, 1, 0, 0, 0, 1, 0, 0, 1])
bwt.insert([0, 1, 1, 1, 0, 1, 1, 0, 1])
bwt.insert([1, 1, 1, 0, 0, 1, 0, 1, 1])
bwt.insert([1, 0, 0, 1, 0, 1, 1, 0, 1])
# # print('nice')
# bwt.print_sorting()
# bwt.print_divergence()

print(bwt.longest_prefix([0, 1, 1, 0, 0, 1, 0, 0, 1]))
print(bwt.longest_prefix([1, 1, 1, 0, 0, 1, 0, 0, 1]))
print(bwt.fastLS([1, 1, 1, 0, 0, 1, 0, 0, 1]))


bwt = DPBWT([1, 2, 3, 4, 5, 6, 7, 8], [1., 2., 3., 4., 5., 6., 7., 8.], 0.5)
bwt.insert([0] * 8)
bwt.insert([1] * 8)
bwt.print_sorting()
print(bwt.fastLS([1, 1, 1, 1, 0, 0, 0, 0]))
