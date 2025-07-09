from threads_arg import Shuttle
import numpy as np

def test_shuttle():
    # G is 2N x M
    G = np.array([[0, 0, 1, 1, 0, 0, 1, 1],
                  [0, 0, 1, 0, 1, 1, 0, 1],
                  [0, 0, 1, 1, 0, 0, 1, 1],
                  [1, 0, 0, 0, 1, 0, 1, 1],
                  [0, 1, 1, 0, 0, 1, 0, 1]])

    targets = [0, 0, 1, 2, 3]
    first_carriers = [3, 4, 0, 0, 1, 1, 0, 0]
    site_hets = [{3, 4}, {4}, {3}, {1, 3, 4}, {2, 4}, {2, 3}, {1}, set()]
    site_positions = [1, 2, 3, 4, 5, 6, 7, 8]

    weft_ids       = [4, 2, 3, 4, 3]
    weft_targets   = [2, 0, 1, 1, 2]
    weft_positions = [3, 4, 5, 6, 7]

    ss = Shuttle(targets, first_carriers, site_hets, site_positions, weft_ids, weft_targets, weft_positions)
    # breakpoint()

    _, ncol = G.shape
    for i in np.arange(ncol):
        expected_carriers = sorted(G[:, i].nonzero()[0].tolist())
        found_carriers = sorted(list(ss.current_carriers))
        assert expected_carriers == found_carriers
        if i < ncol - 1:
            ss.proceed_to_next_site()

def test_GU_mul():
    # G is 2N x M
    G = np.array([[0, 0, 1, 1, 0, 0, 1, 1],
                  [0, 0, 1, 0, 1, 1, 0, 1],
                  [0, 0, 1, 1, 0, 0, 1, 1],
                  [1, 0, 0, 0, 1, 0, 1, 1],
                  [0, 1, 1, 0, 0, 1, 0, 1],
                  [0, 1, 1, 0, 0, 1, 0, 1]])

    targets = [0, 0, 1, 2, 3, 4]
    first_carriers = [3, 4, 0, 0, 1, 1, 0, 0]
    site_hets = [{3, 4}, {4}, {3}, {1, 3, 4}, {2, 4}, {2, 3}, {1}, set()]
    site_positions = [1, 2, 3, 4, 5, 6, 7, 8]

    weft_ids       = [4, 2, 3, 4, 3]
    weft_targets   = [2, 0, 1, 1, 2]
    weft_positions = [3, 4, 5, 6, 7]

    ss = Shuttle(targets, first_carriers, site_hets, site_positions, weft_ids, weft_targets, weft_positions)
    rng = np.random.default_rng(123)
    U = rng.random((8, 4)).astype(np.float32)
    out = np.zeros((3, 4), order="F").astype(np.float32)

    ss.GU_mult(U, out)

    expected_tmp = np.dot(G, U)
    expected = expected_tmp[::2] + expected_tmp[1::2]
    assert np.isclose(expected, out).all()
    # Shuttle(std::vector<int> init_targets,
    #         std::vector<int> _site_first_carriers,
    #         std::vector<std::unordered_set<int>>& _site_hets,
    #         std::vector<int>& _site_positions,
    #         const std::vector<int>& weft_ids,
    #         const std::vector<int>& weft_targets,
    #         const std::vector<int>& weft_positions);
