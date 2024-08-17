import numpy as np
import pytest

import threads_arg

def dummy_threading_instructions():
    return {
        # the "threads" are a list of 3-tuples lists of segment starts, threading targets and segment ages
        "threads": [(np.array([]), np.array([]), np.array([])), 
                    (np.array([4, 8]), np.array([0, 0]), np.array([2.5, 3.0])),
                    (np.array([4, 10, 15]), np.array([0, 1, 0]), np.array([1.0, 3.0, 6.0]))],
        "samples": np.array([0, 1, 2]),
        "positions": np.array([4, 6, 8, 10, 15]),
        "arg_range": np.array([4, 20]),
        "mutations": [np.array([0, 1, 3]),
                      np.array([0, 2, 4]),
                      np.array([1, 3])]
    }

def threading_instructions_are_equal(t1, t2):
    try:
        assert len(t1["threads"]) == len(t2["threads"])

        for (bp1, id1, age1), (bp2, id2, age2) in zip(t1["threads"], t2["threads"]):
            assert (bp1 == bp2).all()
            assert (id1 == id2).all()
            assert (age1 == age2).all()

        for mut1, mut2 in zip(t1["mutations"], t2["mutations"]):
            assert (mut1 == mut2).all()

        assert (t1["samples"] == t2["samples"]).all()
        assert (t1["arg_range"] == t2["arg_range"]).all()
        assert (t1["positions"] == t2["positions"]).all()
        return True
    except AssertionError:
        return False

def test_slice():
    """
    Test the threading instructions slice function
    """
    tt = dummy_threading_instructions()

    # Test illegal intervals
    with pytest.raises(ValueError):
        sliced = threads_arg.utils.slice_threading_instructions(tt, 3, 4)
    with pytest.raises(ValueError):
        sliced = threads_arg.utils.slice_threading_instructions(tt, 4, 4)
    with pytest.raises(ValueError):
        sliced = threads_arg.utils.slice_threading_instructions(tt, 4, 21)
    
    # Test identity on full interval
    # import pdb
    # pdb.set_trace()
    assert threading_instructions_are_equal(threads_arg.utils.slice_threading_instructions(tt, 4, 20), tt)

    # Test interval starting at arg range start and ending between sites
    expected = {
        "threads": [(np.array([]), np.array([]), np.array([])), 
                    (np.array([4]), np.array([0]), np.array([2.5])),
                    (np.array([4]), np.array([0]), np.array([1.0]))],
        "samples": np.array([0, 1, 2]),
        "positions": np.array([4, 6, 8, 10, 15]),
        "arg_range": np.array([4, 7]),
        "mutations": [np.array([0, 1]),
                      np.array([0]),
                      np.array([1])]
    }
    assert threading_instructions_are_equal(threads_arg.utils.slice_threading_instructions(tt, 4, 7), expected)
   
    # Test interval starting at arg range start and ending on a site 
    expected = {
        "threads": [(np.array([]), np.array([]), np.array([])), 
                    (np.array([4, 8]), np.array([0, 0]), np.array([2.5, 3.0])),
                    (np.array([4]), np.array([0]), np.array([1.0]))],
        "samples": np.array([0, 1, 2]),
        "positions": np.array([4, 6, 8, 10, 15]),
        "arg_range": np.array([4, 8]),
        "mutations": [np.array([0, 1]),
                    np.array([0, 2]),
                    np.array([1])]
    }
    assert threading_instructions_are_equal(threads_arg.utils.slice_threading_instructions(tt, 4, 8), expected)

    # Test interval starting at a site and ending at the arg range end
    expected = {
        "threads": [(np.array([]), np.array([]), np.array([])), 
                    (np.array([6, 8]), np.array([0, 0]), np.array([2.5, 3.0])),
                    (np.array([6, 10, 15]), np.array([0, 1, 0]), np.array([1.0, 3.0, 6.0]))],
        "samples": np.array([0, 1, 2]),
        "positions": np.array([4, 6, 8, 10, 15]),
        "arg_range": np.array([6, 20]),
        "mutations": [np.array([1, 3]),
                    np.array([2, 4]),
                    np.array([1, 3])]
    }
    assert threading_instructions_are_equal(threads_arg.utils.slice_threading_instructions(tt, 6, 20), expected)

    # Test interval starting and ending between sites
    expected = {
        "threads": [(np.array([]), np.array([]), np.array([])), 
                    (np.array([7, 8]), np.array([0, 0]), np.array([2.5, 3.0])),
                    (np.array([7, 10, 15]), np.array([0, 1, 0]), np.array([1.0, 3.0, 6.0]))],
        "samples": np.array([0, 1, 2]),
        "positions": np.array([4, 6, 8, 10, 15]),
        "arg_range": np.array([7, 19]),
        "mutations": [np.array([3]),
                    np.array([2, 4]),
                    np.array([3])]
    }
    assert threading_instructions_are_equal(threads_arg.utils.slice_threading_instructions(tt, 7, 19), expected)
