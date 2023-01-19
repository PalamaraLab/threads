from threads import Threads
import numpy as np

def test_mutations():
    mu = 0.001
    rho = 0.001

    # Test 1
    pos1 = np.array([1, 2, 3, 4, 5, 6, 7, 8, 9])
    t1 = Threads(pos1, 100 * rho * pos1, mu, Ne=1)
    t1.insert([0, 1, 0, 0, 0, 1, 0, 0, 1])
    t1.insert([0, 1, 1, 1, 0, 1, 1, 0, 1])
    t1.insert([1, 1, 1, 0, 0, 1, 0, 1, 1])
    t1.insert([1, 0, 0, 1, 0, 1, 1, 0, 1])

    bp_starts_1, target_IDs_1, _, het_sites_1 = t1.thread([1, 1, 1, 0, 0, 1, 0, 1, 0])
    assert bp_starts_1 == [1]
    assert target_IDs_1 == [2]
    assert het_sites_1 == [9]
    
    # Test 2
    pos2 = np.array([1, 2, 3, 4, 5, 6, 7, 8])
    t2 = Threads(pos2, 100 * rho * pos2, mu, Ne=1)
    t2.insert([0] * 8)
    t2.insert([1] * 8)
    bp_starts_2, target_IDs_2, _, het_sites_2 = t2.thread([1, 1, 0, 0, 0, 0, 1, 1])
    assert bp_starts_2 == [1, 3, 7]
    assert target_IDs_2 == [1, 0, 1]
    assert het_sites_2 == []

    # Test 3
    pos3 = np.array([1, 2, 3, 4, 5, 6, 7, 8])
    t3 = Threads(pos3, 100 * rho * pos3, mu, Ne=1)
    t3.insert([0, 0, 0, 0, 1, 0, 0, 0])
    t3.insert([0, 0, 0, 0, 1, 0, 0, 0])
    t3.insert([1, 1, 1, 1, 0, 0, 0, 0])
    t3.insert([1, 1, 1, 1, 0, 0, 0, 0])
    bp_starts_3, target_IDs_3, _, het_sites_3 = t3.thread([1, 1, 1, 1, 1, 0, 0, 0])
    assert bp_starts_3 == [1]
    assert target_IDs_3 in [[2], [3]]
    assert het_sites_3 == [5]
