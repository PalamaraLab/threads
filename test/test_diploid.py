from threads import Threads
import numpy as np

# mu = 0.001
# rho = 0.001
# # Test 1: perfect copying
# pos1 = np.array([1, 2, 3, 4, 5, 6, 7, 8, 9])
# t1 = Threads(pos1, 100 * rho * pos1, mu, Ne=1)
# G1 = np.array([[0, 1, 0, 0, 0, 1, 0, 0, 1],
#                 [0, 1, 1, 1, 0, 1, 1, 0, 1],
#                 [1, 1, 1, 0, 0, 1, 0, 1, 1],
#                 [1, 0, 0, 1, 0, 1, 1, 0, 1]])
# t1.insert(G1[0])
# t1.insert(G1[1])
# t1.insert(G1[2])
# t1.insert(G1[3])
# path_a, path_b = t1.diploid_ls(G1[1] + G1[2])
# phased = t1.phase(G1[1] + G1[2])

def test_diploid():
    mu = 0.001
    rho = 0.001
    # Test 1: perfect copying
    pos1 = np.array([1, 2, 3, 4, 5, 6, 7, 8, 9])
    t1 = Threads(pos1, 100 * rho * pos1, mu, Ne=1)
    G1 = np.array([[0, 1, 0, 0, 0, 1, 0, 0, 1],
                   [0, 1, 1, 1, 0, 1, 1, 0, 1],
                   [1, 1, 1, 0, 0, 1, 0, 1, 1],
                   [1, 0, 0, 1, 0, 1, 1, 0, 1]])
    t1.insert(G1[0])
    t1.insert(G1[1])
    t1.insert(G1[2])
    t1.insert(G1[3])
    path_a, path_b = t1.diploid_ls(G1[1] + G1[2])
    import pdb
    pdb.set_trace()
    assert len(path_a) == len(path_b) == 1
    assert path_a[0][1][0] == 1 and path_b[0][1][0] == 2 or path_a[0][1][0] == 2 and path_b[0][1][0] == 1

    # todo more tests
    # Test 2: copying with one error
    pos2 = np.array([1, 2, 3, 4, 5, 6, 7, 8, 9])
    t2 = Threads(pos1, 100 * rho * pos2, mu, Ne=1)
    G2 = np.array([[0, 1, 0, 0, 0, 1, 0, 0, 1],
                [0, 1, 1, 1, 0, 1, 1, 0, 1],
                [1, 1, 1, 0, 0, 1, 0, 1, 1],
                [1, 0, 0, 1, 0, 1, 1, 0, 1]])
    t2.insert(G2[0])
    t2.insert(G2[1])
    t2.insert(G2[2])
    t2.insert(G2[3])
    hap_a = G2[0].copy()
    hap_a[2] = 1
    hap_b = G2[1]

    path_a, path_b = t2.fastLS_diploid(hap_a + hap_b)
    assert len(path_a) == len(path_b) == 1
    assert path_a[0][1][0] == 1 and path_b[0][1][0] == 0 or path_a[0][1][0] == 0 and path_b[0][1][0] == 1

    # Test 3: one recombination
    pos3 = np.array([1, 2, 3, 4, 5, 6, 7, 8, 9])
    t3 = Threads(pos1, 10 * pos3, mu, Ne=1)
    G3 = np.array([[0, 1, 0, 0, 1, 1, 0, 0, 1],
                [0, 1, 1, 1, 1, 1, 1, 0, 1],
                [1, 1, 1, 0, 0, 1, 0, 1, 1],
                [1, 0, 0, 1, 0, 0, 1, 1, 1]])
    t3.insert(G3[0])
    t3.insert(G3[1])
    t3.insert(G3[2])
    t3.insert(G3[3])
    hap_a = G3[0].copy()
    hap_a[4:] = G3[2, 4:]
    hap_b = G3[1]
    path_a, path_b = t3.fastLS_diploid(hap_a + hap_b)
    assert len(path_a) == 1 and len(path_b) == 2 or len(path_a) == 2 and len(path_b) == 1
    if len(path_a) == 2:
        assert path_a[0][1][0] == 0 
        assert path_a[1][1][0] == 2
        assert path_b[0][1][0] == 1
    else:
        assert path_b[0][1][0] == 0 
        assert path_b[1][1][0] == 2
        assert path_a[0][1][0] == 1
