"""
Tests for normalization.py: Normalizer class (demography-based TMRCA correction).
Uses a minimal synthetic demography and small sample count to keep tests fast.
"""
import tempfile

import numpy as np
import pytest

from threads_arg import ThreadingInstructions


# ---------------------------------------------------------------------------
# Fixtures
# ---------------------------------------------------------------------------
@pytest.fixture(scope="module")
def demo_file():
    """Minimal constant-Ne demography file."""
    with tempfile.NamedTemporaryFile(mode="w", suffix=".demo", delete=False) as f:
        f.write("0.0\t10000\n")
        f.write("500.0\t10000\n")
        f.flush()
        return f.name


@pytest.fixture(scope="module")
def normalizer(demo_file):
    from threads_arg.normalization import Normalizer
    return Normalizer(demo_file, num_samples=20)


@pytest.fixture(scope="module")
def simple_instructions():
    """Small ThreadingInstructions: 20 haploids, 5 sites."""
    # Sample 0 is the reference (no threading info needed for it,
    # but ThreadingInstructions expects entries for all samples)
    starts = [[0]] * 20
    tmrcas = [[float(i * 50 + 100)] for i in range(20)]  # 100, 150, ..., 1050
    targets = [[0]] * 20
    mismatches = [[]] * 20
    positions = [100, 200, 300, 400, 500]
    return ThreadingInstructions(starts, tmrcas, targets, mismatches, positions, 0, 600)


# ===================================================================
# Normalizer construction
# ===================================================================
class TestNormalizerInit:
    def test_num_samples_stored(self, normalizer):
        assert normalizer.num_samples == 20

    def test_demography_created(self, normalizer):
        assert normalizer.demography is not None


# ===================================================================
# Normalizer.simulation
# ===================================================================
class TestSimulation:
    def test_returns_tree_sequence(self, normalizer):
        ts = normalizer.simulation(1e6, random_seed=42)
        assert ts.num_samples == 20  # 10 diploid = 20 haploid

    def test_has_internal_nodes(self, normalizer):
        ts = normalizer.simulation(1e6, random_seed=42)
        internal_times = [n.time for n in ts.nodes() if n.time > 0]
        assert len(internal_times) > 0

    def test_deterministic_with_seed(self, normalizer):
        ts1 = normalizer.simulation(1e6, random_seed=99)
        ts2 = normalizer.simulation(1e6, random_seed=99)
        times1 = sorted([n.time for n in ts1.nodes()])
        times2 = sorted([n.time for n in ts2.nodes()])
        assert times1 == times2


# ===================================================================
# Normalizer.normalize
# ===================================================================
class TestNormalize:
    def test_preserves_starts(self, normalizer, simple_instructions):
        result = normalizer.normalize(simple_instructions, num_seeds=10)
        assert result.all_starts() == simple_instructions.all_starts()

    def test_preserves_targets(self, normalizer, simple_instructions):
        result = normalizer.normalize(simple_instructions, num_seeds=10)
        assert result.all_targets() == simple_instructions.all_targets()

    def test_preserves_mismatches(self, normalizer, simple_instructions):
        result = normalizer.normalize(simple_instructions, num_seeds=10)
        assert result.all_mismatches() == simple_instructions.all_mismatches()

    def test_preserves_positions(self, normalizer, simple_instructions):
        result = normalizer.normalize(simple_instructions, num_seeds=10)
        assert result.positions == simple_instructions.positions

    def test_preserves_num_samples(self, normalizer, simple_instructions):
        result = normalizer.normalize(simple_instructions, num_seeds=10)
        assert result.num_samples == simple_instructions.num_samples

    def test_preserves_region(self, normalizer, simple_instructions):
        result = normalizer.normalize(simple_instructions, num_seeds=10)
        assert result.start == simple_instructions.start
        assert result.end == simple_instructions.end

    def test_tmrcas_changed(self, normalizer, simple_instructions):
        """Normalized TMRCAs should differ from originals (unless by coincidence)."""
        result = normalizer.normalize(simple_instructions, num_seeds=10)
        old = [h for hs in simple_instructions.all_tmrcas() for h in hs]
        new = [h for hs in result.all_tmrcas() for h in hs]
        assert len(old) == len(new)
        # At least some should differ
        n_changed = sum(1 for a, b in zip(old, new) if abs(a - b) > 1e-6)
        assert n_changed > 0

    def test_tmrcas_non_negative(self, normalizer, simple_instructions):
        result = normalizer.normalize(simple_instructions, num_seeds=10)
        for tmrca_vec in result.all_tmrcas():
            for t in tmrca_vec:
                assert t >= 0, f"negative TMRCA: {t}"

    def test_monotonic_mapping(self, normalizer, simple_instructions):
        """If original TMRCA_a < TMRCA_b, normalized should preserve order."""
        result = normalizer.normalize(simple_instructions, num_seeds=10)
        old_tmrcas = [hs[0] for hs in simple_instructions.all_tmrcas()]
        new_tmrcas = [hs[0] for hs in result.all_tmrcas()]
        # Sort by old, check new is non-decreasing
        pairs = sorted(zip(old_tmrcas, new_tmrcas))
        new_sorted = [p[1] for p in pairs]
        for i in range(1, len(new_sorted)):
            assert new_sorted[i] >= new_sorted[i - 1] - 1e-10, \
                f"monotonicity violated: {new_sorted[i-1]} > {new_sorted[i]}"

    def test_sample_count_mismatch_raises(self, demo_file):
        """Normalizer should assert if num_samples doesn't match instructions."""
        from threads_arg.normalization import Normalizer
        norm = Normalizer(demo_file, num_samples=10)
        # Instructions with 20 samples but normalizer expects 10
        starts = [[0]] * 20
        tmrcas = [[100.0]] * 20
        targets = [[0]] * 20
        mismatches = [[]] * 20
        ti = ThreadingInstructions(starts, tmrcas, targets, mismatches, [100], 0, 200)
        with pytest.raises(AssertionError):
            norm.normalize(ti, num_seeds=5)
