from threads import Threads
import numpy as np

def test_trimmed_positions():
    # Test 1
    pos1 = np.array([1, 2, 3, 4, 5, 6, 7, 8, 9])
    t1 = Threads(pos1, pos1, 1., Ne=1, burn_in_left=3, burn_in_right=1)
    assert t1.threading_start == 4.
    assert t1.threading_end == 8.
    assert t1.trimmed_positions() == [4., 5., 6., 7., 8.]


    pos2 = np.array([1, 2, 4, 6, 8, 10, 12, 14])
    t2 = Threads(pos2, pos2, 1., Ne=1, burn_in_left=2, burn_in_right=1)
    assert t2.threading_start == 3.
    assert t2.threading_end == 13.
    assert t2.trimmed_positions() == [4, 6, 8, 10, 12]

def test_burn_in():
    mu = 0.001
    rho = 0.001

    pos1 = np.array([1, 2, 3, 4, 5, 6, 7, 8, 9])
    t1 = Threads(pos1, 100 * rho * pos1, mu, Ne=1, burn_in_left=2, burn_in_right=1)

    t1.insert([0, 1, 0, 0, 0, 1, 0, 0, 1])
    t1.insert([0, 1, 1, 1, 0, 1, 1, 0, 1])
    t1.insert([1, 1, 1, 0, 0, 1, 0, 1, 1])
    t1.insert([1, 0, 0, 1, 0, 1, 1, 0, 1])

    bp_starts_1, target_IDs_1, _, het_sites_1 = t1.thread([1, 1, 1, 0, 0, 1, 0, 1, 0])
    assert bp_starts_1 == [3]
    assert target_IDs_1 == [2]
    assert het_sites_1 == []
    
    # Test 2
    pos2 = np.array([0, 2, 4, 6, 8, 10, 12, 14])
    t2 = Threads(pos2, 100 * rho * pos2, mu, Ne=1, burn_in_left=3, burn_in_right=1)
    t2.insert([0] * 8)
    t2.insert([1] * 8)
    bp_starts_2, target_IDs_2, _, het_sites_2 = t2.thread([1, 1, 0, 0, 0, 0, 1, 1])
    assert bp_starts_2 == [3, 11]
    assert target_IDs_2 == [0, 1]
    assert het_sites_2 == []

    # Test 3
    pos3 = np.array([1, 2, 3, 4, 5, 6, 7, 8])
    t3 = Threads(pos3, 100 * rho * pos3, mu, Ne=1, burn_in_left=0, burn_in_right=4)
    t3.insert([0, 0, 0, 0, 1, 0, 0, 0])
    t3.insert([0, 0, 0, 0, 1, 0, 0, 0])
    t3.insert([1, 1, 1, 1, 0, 0, 0, 0])
    t3.insert([1, 1, 1, 1, 0, 0, 0, 0])
    bp_starts_3, target_IDs_3, _, het_sites_3 = t3.thread([1, 1, 1, 1, 1, 0, 0, 0])
    assert bp_starts_3 == [1]
    assert target_IDs_3 in [[2], [3]]
    assert het_sites_3 == []
