import numpy as np
from threads import TGEN, Threads

def test_small_tgen():
    mu = 0.001
    rho = 0.001

    G = np.array([
        [0, 1, 0, 0, 0, 1, 0, 0, 1],
        [0, 1, 1, 1, 0, 1, 1, 0, 1],
        [1, 1, 1, 0, 0, 1, 0, 1, 1],
        [1, 0, 0, 1, 0, 1, 1, 0, 1],
        [1, 1, 1, 0, 0, 1, 0, 1, 0]
    ], dtype=int)

    # Test 1
    pos = np.array([1, 2, 3, 4, 5, 6, 7, 8, 9])
    t = Threads(pos, 100 * rho * pos, mu, Ne=1)
    t.insert(G[0])
    bp_starts_1, target_ids_1, _, het_sites_1 = [], [[]], [], pos[np.nonzero(G[0])[0].tolist()]
    bp_starts_2, target_ids_2, _, het_sites_2 = t.thread(G[1])
    bp_starts_3, target_ids_3, _, het_sites_3 = t.thread(G[2])
    bp_starts_4, target_ids_4, _, het_sites_4 = t.thread(G[3])
    bp_starts_5, target_ids_5, _, het_sites_5 = t.thread(G[4])

    bp_starts = [bp_starts_1, bp_starts_2, bp_starts_3, bp_starts_4, bp_starts_5]
    print( [target_ids_1, target_ids_2, target_ids_3, target_ids_4, target_ids_5])
    target_ids = [[]] + [[t[0] for t in t_list] for t_list in [target_ids_2, target_ids_3, target_ids_4, target_ids_5]]
    het_sites = [het_sites_1, het_sites_2, het_sites_3, het_sites_4, het_sites_5]

    tgen = TGEN(pos, bp_starts, target_ids, het_sites)
    tgen.query(1, 3, [2, 3])
    assert (tgen.query(1, 3, [2, 3]) == G[[2, 3], 0:3]).all()
    assert (tgen.query(5, 9, [1, 2, 4]) == G[[1, 2, 4], 4:]).all()
    assert (tgen.query(1, 9, [0, 1, 2, 3, 4]) == G).all()

def test_big_tgen():
    N = 300
    M = 1000
    rng = np.random.default_rng(130222)
    big_matrix = rng.random((N, M)) - 0.5
    big_matrix[big_matrix < 0] = 0
    big_matrix[big_matrix > 0] = 1
    big_matrix = big_matrix.astype(int)

    mu = 0.001
    rho = 0.001
    pos = np.arange(M)
    t = Threads(pos, 100 * rho * pos, mu, Ne=1)

    bp_starts, target_ids, het_sites = [], [], []
    for i in range(N):
        if i == 0:
            t.insert(big_matrix[0])
            bp_starts.append([])
            target_ids.append([])
            het_sites.append(pos[np.nonzero(big_matrix[0])[0].tolist()])
        else:
            bp, targets, _, het = t.thread(big_matrix[i])
            bp_starts.append(bp)
            target_ids.append([t[0] for t in targets])
            het_sites.append(het)
    
    tgen = TGEN(pos, bp_starts, target_ids, het_sites)
    assert (tgen.query(pos[0], pos[-1], np.arange(N)) == big_matrix).all()
    assert (tgen.query(2, 5, [12, 0, 2]) == big_matrix[[12, 0, 2], 2:6]).all()
