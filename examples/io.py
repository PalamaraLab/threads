from TDPBWT import DPBWT
import numpy as np
import gzip
mu = 0.001
rho = 0.001

# Test 1
pos1 = np.array([1, 2, 3, 4, 5, 6, 7, 8, 9])
bwt = DPBWT(pos1, 100 * rho * pos1, mu, Ne=1)
haps = np.array([[0, 1, 0, 0, 0, 1, 0, 0, 1],
                 [0, 1, 1, 1, 0, 1, 1, 0, 1],
                 [1, 1, 1, 0, 0, 1, 0, 1, 1],
                 [1, 0, 0, 1, 0, 1, 1, 0, 1],
                 [1, 1, 1, 0, 0, 1, 0, 0, 1]])

with gzip.open("example_data/hack.gz", "w+") as hack_gz:
    for hap in haps:
        line_str = " ".join([str(j) for j in hap]) + "\n"
        hack_gz.write(line_str.encode())

threads = bwt.thread_from_file("example_data/hack.gz", n_cycle=3)
print(threads)
