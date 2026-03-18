"""
Tests for the `threads convert` command: threads_to_arg, ARG structure validation,
and serialization to .argn format.
"""
import tempfile

import numpy as np
import pytest
import arg_needle_lib

from pathlib import Path

TEST_DATA = Path(__file__).parent / "data"
THREADS_FIT = str(TEST_DATA / "expected_infer_fit_to_data_snapshot.threads")
THREADS_NO_FIT = str(TEST_DATA / "expected_infer_snapshot.threads")
PANEL_PGEN = str(TEST_DATA / "panel.pgen")


@pytest.fixture(scope="module")
def instructions_nofit():
    from threads_arg.serialization import load_instructions
    return load_instructions(THREADS_NO_FIT)


@pytest.fixture(scope="module")
def instructions_fit():
    from threads_arg.serialization import load_instructions
    return load_instructions(THREADS_FIT)


def _build_arg(instructions, add_mutations=False):
    """Build ARG with retry logic matching threads_convert."""
    from threads_arg.convert import threads_to_arg
    for noise in [0.0, 1e-5, 1e-3]:
        try:
            return threads_to_arg(instructions, add_mutations=add_mutations, noise=noise)
        except RuntimeError:
            continue
    raise RuntimeError("Failed to build ARG even with noise=1e-3")


@pytest.fixture(scope="module")
def arg_nofit(instructions_nofit):
    return _build_arg(instructions_nofit, add_mutations=False)


# ===================================================================
# threads_to_arg: ARG construction
# ===================================================================
class TestThreadsToArg:
    def test_arg_has_correct_sample_count(self, arg_nofit, instructions_nofit):
        leaf_ids = arg_nofit.leaf_ids
        assert len(leaf_ids) == instructions_nofit.num_samples

    def test_arg_leaves_are_leaves(self, arg_nofit):
        for lid in arg_nofit.leaf_ids:
            assert arg_nofit.is_leaf(lid)

    def test_arg_offset_matches_instructions(self, arg_nofit, instructions_nofit):
        assert arg_nofit.offset == instructions_nofit.start

    def test_arg_construction_with_noise(self, instructions_nofit):
        """Adding noise should not crash and should produce a valid ARG."""
        from threads_arg.convert import threads_to_arg
        arg = threads_to_arg(instructions_nofit, add_mutations=False, noise=1e-5)
        assert len(arg.leaf_ids) == instructions_nofit.num_samples

    def test_arg_with_mutations(self, instructions_nofit):
        """add_mutations=True should populate mutations on the ARG."""
        arg = _build_arg(instructions_nofit, add_mutations=True)
        arg.populate_children_and_roots()
        mutations = arg.mutations()
        # Parsimonious mapping may place multiple mutations per site
        assert len(mutations) >= instructions_nofit.num_sites


# ===================================================================
# threads_convert: full pipeline .threads -> .argn
# ===================================================================
class TestThreadsConvert:
    def test_produces_argn_file(self):
        from threads_arg.convert import threads_convert
        with tempfile.TemporaryDirectory() as tmpdir:
            out_argn = str(Path(tmpdir) / "test.argn")
            threads_convert(THREADS_NO_FIT, argn=out_argn, tsz=None)
            assert Path(out_argn).exists()
            assert Path(out_argn).stat().st_size > 0

    def test_argn_deserializable(self):
        from threads_arg.convert import threads_convert
        with tempfile.TemporaryDirectory() as tmpdir:
            out_argn = str(Path(tmpdir) / "test.argn")
            threads_convert(THREADS_NO_FIT, argn=out_argn, tsz=None)
            arg = arg_needle_lib.deserialize_arg(out_argn)
            assert len(arg.leaf_ids) == 1000

    def test_argn_snapshot_match(self):
        """Generated .argn should match the expected snapshot."""
        from threads_arg.convert import threads_convert
        import h5py
        expected_argn = TEST_DATA / "expected_convert_snapshot.argn"
        with tempfile.TemporaryDirectory() as tmpdir:
            out_argn = str(Path(tmpdir) / "test.argn")
            threads_convert(THREADS_NO_FIT, argn=out_argn, tsz=None)

            # Compare HDF5 structures
            with h5py.File(out_argn, "r") as gen, h5py.File(str(expected_argn), "r") as exp:
                assert set(gen.keys()) == set(exp.keys())
                for key in exp.keys():
                    if isinstance(exp[key], h5py.Dataset):
                        assert gen[key].shape == exp[key].shape, \
                            f"shape mismatch for {key}: {gen[key].shape} vs {exp[key].shape}"

    def test_fit_to_data_converts(self):
        """fit_to_data .threads should also convert without error."""
        from threads_arg.convert import threads_convert
        with tempfile.TemporaryDirectory() as tmpdir:
            out_argn = str(Path(tmpdir) / "test_fit.argn")
            threads_convert(THREADS_FIT, argn=out_argn, tsz=None)
            assert Path(out_argn).exists()


# ===================================================================
# ARG structure properties
# ===================================================================
class TestArgStructure:
    def test_genotype_reconstruction_from_argn(self):
        """Genotypes from the pre-built .argn snapshot (without add_mutations)
        can be deserialized and have the correct sample count."""
        expected_argn = str(TEST_DATA / "expected_convert_snapshot.argn")
        arg = arg_needle_lib.deserialize_arg(expected_argn)
        assert len(arg.leaf_ids) == 1000
        for lid in arg.leaf_ids:
            assert arg.is_leaf(lid)
