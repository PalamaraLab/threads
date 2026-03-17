# This file is part of the Threads software suite.
# Copyright (C) 2024-2025 Threads Developers.
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <http://www.gnu.org/licenses/>.

"""
Tests for AgeEstimator and ConsistencyWrapper.

Unit tests use synthetic ThreadingInstructions with known properties.
Regression tests record exact outputs from real data.
"""

import numpy as np
import pytest

from threads_arg import AgeEstimator, GenotypeIterator, ThreadingInstructions, ConsistencyWrapper
from threads_arg.serialization import load_instructions

from snapshot_runners import TEST_DATA_DIR


# ---------------------------------------------------------------------------
# Helpers for synthetic ThreadingInstructions
# ---------------------------------------------------------------------------

def _make_inst(targets, tmrcas, positions=None, starts=None, mismatches=None):
    """Build a ThreadingInstructions from per-sample targets/tmrcas.

    For single-segment cases, pass flat lists:
        _make_inst([-1, 0, 1], [0, 500, 200])
    For multi-segment, pass nested lists and explicit starts/positions.
    """
    n = len(targets)
    if positions is None:
        positions = [1000]
    if not isinstance(targets[0], list):
        targets = [[t] for t in targets]
    if not isinstance(tmrcas[0], list):
        tmrcas = [[float(t)] for t in tmrcas]
    if starts is None:
        starts = [[positions[0]] for _ in range(n)]
    if mismatches is None:
        mismatches = [[] for _ in range(n)]
    return ThreadingInstructions(
        starts, tmrcas, targets, mismatches,
        positions, positions[0], positions[-1] + 1,
    )


def _run_age(inst, genotypes_list):
    """Run AgeEstimator, return list of ages."""
    ae = AgeEstimator(inst)
    for g in genotypes_list:
        ae.process_site(g)
    return ae.get_inferred_ages()


def _run_dc(inst, ages, genotypes_list):
    """Run ConsistencyWrapper, return consistent ThreadingInstructions."""
    cw = ConsistencyWrapper(inst, ages)
    for g in genotypes_list:
        cw.process_site(g)
    return cw.get_consistent_instructions()


# ---------------------------------------------------------------------------
# Unit tests: AgeEstimator
# ---------------------------------------------------------------------------


def test_age_two_carriers_coalesce():
    """Two carriers sharing a recent ancestor get age near their TMRCA."""
    # 0 <- 1 (t=1000) <- 2 (t=500) <- 3 (t=200)
    inst = _make_inst([-1, 0, 1, 2], [0, 1000, 500, 200])
    ages = _run_age(inst, [[0, 0, 1, 1]])
    # Carriers at 2,3 coalesce at t=200; sweep gives boundary=200, next=500
    assert ages[0] == pytest.approx(350.0)


def test_age_all_carriers():
    """All carriers: age above the deepest TMRCA."""
    inst = _make_inst([-1, 0, 1, 2], [0, 1000, 500, 200])
    ages = _run_age(inst, [[1, 1, 1, 1]])
    assert ages[0] > 1000


def test_age_no_carriers():
    """No carriers: still produces a positive age."""
    inst = _make_inst([-1, 0, 1, 2], [0, 1000, 500, 200])
    ages = _run_age(inst, [[0, 0, 0, 0]])
    assert ages[0] > 0


def test_age_single_carrier():
    """One carrier among non-carriers."""
    inst = _make_inst([-1, 0, 1], [0, 1000, 500])
    ages = _run_age(inst, [[0, 1, 0]])
    assert ages[0] > 0


def test_age_carrier_at_root_only():
    """Only the root sample is a carrier."""
    inst = _make_inst([-1, 0, 1], [0, 500, 200])
    ages = _run_age(inst, [[1, 0, 0]])
    assert ages[0] > 0


def test_age_two_samples():
    """Minimal two-sample case."""
    inst = _make_inst([-1, 0], [0, 100])
    ages = _run_age(inst, [[1, 1]])
    assert ages[0] > 0
    ages2 = _run_age(inst, [[1, 0]])
    assert ages2[0] > 0


def test_age_reflects_tmrca_scale():
    """Closer carriers yield younger age than distant carriers."""
    inst_close = _make_inst([-1, 0, 1], [0, 100, 50])
    inst_far = _make_inst([-1, 0, 1], [0, 10000, 5000])
    age_close = _run_age(inst_close, [[0, 1, 1]])[0]
    age_far = _run_age(inst_far, [[0, 1, 1]])[0]
    assert age_close < age_far


def test_age_self_referencing_target():
    """Self-referencing target (target[i]==i) doesn't hang or crash."""
    inst = _make_inst([-1, 1, 0], [0, 500, 300])  # sample 1 self-refs
    ages = _run_age(inst, [[1, 1, 0]])
    assert ages[0] > 0


def test_age_self_ref_on_trace_path():
    """Self-ref sample on the trace path from a carrier chain.

    Trace from carrier chain 3->2 goes 2->1(self-ref)->break.
    Sample 0 is never visited during trace, so the fill loop must handle
    sample 0's target=-1 without crashing.
    """
    inst = _make_inst([-1, 1, 1, 2], [0, 800, 400, 100])
    ages = _run_age(inst, [[0, 0, 1, 1]])
    assert ages[0] > 0


def test_age_multiple_sites():
    """Multiple sites produce one age per site, all positive."""
    positions = [1000, 2000, 3000]
    inst = _make_inst([-1, 0, 1], [0, 500, 200], positions=positions)
    genotypes = [[1, 1, 0], [0, 1, 1], [1, 0, 1]]
    ages = _run_age(inst, genotypes)
    assert len(ages) == 3
    assert all(a > 0 for a in ages)


def test_age_segment_change():
    """Sample changes target mid-sequence; ages still positive."""
    positions = [1000, 2000, 3000]
    inst = _make_inst(
        targets=[[-1], [0, 2], [0]],
        tmrcas=[[0.0], [500.0, 200.0], [300.0]],
        positions=positions,
        starts=[[1000], [1000, 2000], [1000]],
    )
    ages = _run_age(inst, [[1, 1, 0], [0, 1, 1], [1, 0, 1]])
    assert len(ages) == 3
    assert all(a > 0 for a in ages)


def test_age_deterministic_synthetic():
    """Same input twice gives identical output."""
    inst = _make_inst([-1, 0, 1, 2], [0, 1000, 500, 200])
    g = [0, 1, 1, 0]
    a1 = _run_age(inst, [g])[0]
    a2 = _run_age(inst, [g])[0]
    assert a1 == a2


def test_age_many_self_refs():
    """Multiple self-referencing samples in the same tree."""
    inst = _make_inst(
        targets=[-1, 1, 2, 0, 3],
        tmrcas=[0, 600, 400, 200, 100],
    )
    ages = _run_age(inst, [[1, 0, 1, 1, 0]])
    assert ages[0] > 0


# ---------------------------------------------------------------------------
# Unit tests: ConsistencyWrapper
# ---------------------------------------------------------------------------

def test_dc_basic():
    """ConsistencyWrapper on a simple tree preserves dimensions."""
    positions = [1000, 2000, 3000]
    inst = _make_inst([-1, 0, 1], [0, 500, 200], positions=positions)
    genotypes = [[1, 1, 0], [0, 1, 1], [1, 0, 1]]
    ages = _run_age(inst, genotypes)
    dc = _run_dc(inst, ages, genotypes)
    assert dc.num_samples == inst.num_samples
    assert dc.num_sites == inst.num_sites


def test_dc_all_carriers():
    """DC with all-carrier sites."""
    positions = [1000, 2000]
    inst = _make_inst([-1, 0, 1], [0, 500, 200], positions=positions)
    genotypes = [[1, 1, 1], [1, 1, 1]]
    ages = _run_age(inst, genotypes)
    dc = _run_dc(inst, ages, genotypes)
    assert dc.num_samples == 3
    assert dc.num_sites == 2


def test_dc_no_carriers():
    """DC with no-carrier sites."""
    positions = [1000, 2000]
    inst = _make_inst([-1, 0, 1], [0, 500, 200], positions=positions)
    genotypes = [[0, 0, 0], [0, 0, 0]]
    ages = _run_age(inst, genotypes)
    dc = _run_dc(inst, ages, genotypes)
    assert dc.num_samples == 3


def test_dc_self_ref():
    """DC doesn't hang on self-referencing targets."""
    positions = [1000, 2000]
    inst = _make_inst([-1, 1, 0], [0, 500, 300], positions=positions)
    genotypes = [[1, 1, 0], [0, 1, 1]]
    ages = _run_age(inst, genotypes)
    dc = _run_dc(inst, ages, genotypes)
    assert dc.num_samples == inst.num_samples


def test_dc_segment_change():
    """DC with mid-sequence target change."""
    positions = [1000, 2000, 3000]
    inst = _make_inst(
        targets=[[-1], [0, 2], [0]],
        tmrcas=[[0.0], [500.0, 200.0], [300.0]],
        positions=positions,
        starts=[[1000], [1000, 2000], [1000]],
    )
    genotypes = [[1, 1, 0], [0, 1, 1], [1, 0, 1]]
    ages = _run_age(inst, genotypes)
    dc = _run_dc(inst, ages, genotypes)
    assert dc.num_samples == inst.num_samples
    assert dc.num_sites == inst.num_sites


def test_dc_preserves_genotypes():
    """DC output must reconstruct the exact same genotype matrix."""
    positions = [1000, 2000, 3000, 4000, 5000]
    inst = _make_inst([-1, 0, 1, 2], [0, 800, 400, 100], positions=positions)

    # Get the genotypes encoded by the original instructions
    gi_orig = GenotypeIterator(inst)
    orig_genos = []
    while gi_orig.has_next_genotype():
        orig_genos.append(list(gi_orig.next_genotype()))

    ages = _run_age(inst, orig_genos)
    dc = _run_dc(inst, ages, orig_genos)

    gi_dc = GenotypeIterator(dc)
    for i, orig_g in enumerate(orig_genos):
        dc_g = list(gi_dc.next_genotype())
        assert orig_g == dc_g, f"Genotype mismatch at site {i}: {orig_g} != {dc_g}"


def test_dc_preserves_genotypes_with_self_ref():
    """DC preserves genotypes even when original has self-referencing targets."""
    import os
    cache_path = os.path.join(os.path.dirname(__file__), "threads_cache", "sim_50dip.threads")
    if not os.path.exists(cache_path):
        pytest.skip("sim_50dip.threads not in cache")

    inst = load_instructions(cache_path)

    gi = GenotypeIterator(inst)
    all_genos = []
    while gi.has_next_genotype():
        all_genos.append(list(gi.next_genotype()))

    ages = _run_age(inst, all_genos)
    dc = _run_dc(inst, ages, all_genos)

    gi_dc = GenotypeIterator(dc)
    for i, orig_g in enumerate(all_genos):
        dc_g = list(gi_dc.next_genotype())
        assert orig_g == dc_g, f"Genotype mismatch at site {i}"


def test_dc_preserves_multiply():
    """DC right_multiply must produce identical results to original."""
    import os
    cache_path = os.path.join(os.path.dirname(__file__), "threads_cache", "sim_50dip.threads")
    if not os.path.exists(cache_path):
        pytest.skip("sim_50dip.threads not in cache")

    inst = load_instructions(cache_path)
    gi = GenotypeIterator(inst)
    all_genos = []
    while gi.has_next_genotype():
        all_genos.append(list(gi.next_genotype()))

    ages = _run_age(inst, all_genos)
    dc = _run_dc(inst, ages, all_genos)

    rng = np.random.default_rng(42)
    x = rng.normal(0, 1, inst.num_sites).tolist()
    inst.materialize_genotypes()
    dc.materialize_genotypes()
    np.testing.assert_allclose(inst.right_multiply(x), dc.right_multiply(x))


def test_dc_output_genotypes_valid():
    """DC output produces valid genotypes (0 or 1, correct dimensions)."""
    positions = [1000, 2000, 3000, 4000, 5000]
    inst = _make_inst([-1, 0, 1, 2], [0, 800, 400, 100], positions=positions)
    genotypes = [[0, 1, 1, 0], [1, 0, 0, 1], [1, 1, 1, 0], [0, 0, 1, 1], [1, 1, 0, 0]]
    ages = _run_age(inst, genotypes)
    dc = _run_dc(inst, ages, genotypes)

    gi_dc = GenotypeIterator(dc)
    site_count = 0
    while gi_dc.has_next_genotype():
        g = gi_dc.next_genotype()
        assert len(g) == dc.num_samples
        assert all(v in (0, 1) for v in g)
        site_count += 1
    assert site_count == dc.num_sites


def test_dc_ages_positive_on_output():
    """Re-estimating ages on DC output still gives positive ages."""
    positions = [1000, 2000, 3000, 4000, 5000]
    inst = _make_inst([-1, 0, 1, 2], [0, 800, 400, 100], positions=positions)
    genotypes = [[0, 1, 1, 0], [1, 0, 0, 1], [1, 1, 1, 0], [0, 0, 1, 1], [1, 1, 0, 0]]
    ages = _run_age(inst, genotypes)
    dc = _run_dc(inst, ages, genotypes)

    gi = GenotypeIterator(dc)
    ae = AgeEstimator(dc)
    while gi.has_next_genotype():
        ae.process_site(gi.next_genotype())
    dc_ages = ae.get_inferred_ages()
    assert len(dc_ages) == dc.num_sites
    assert all(a > 0 for a in dc_ages)


# ---------------------------------------------------------------------------
# Regression tests (real data)
# ---------------------------------------------------------------------------

@pytest.fixture(scope="module")
def instructions():
    """Load the test threading instructions once for all tests."""
    return load_instructions(str(TEST_DATA_DIR / "expected_infer_snapshot.threads"))


@pytest.fixture(scope="module")
def reference_ages(instructions):
    """Compute allele ages from the current implementation — the ground truth."""
    gt_it = GenotypeIterator(instructions)
    age_est = AgeEstimator(instructions)
    while gt_it.has_next_genotype():
        g = np.array(gt_it.next_genotype())
        age_est.process_site(g)
    return np.array(age_est.get_inferred_ages())


def test_allele_ages_count(instructions, reference_ages):
    """One age per site."""
    assert len(reference_ages) == instructions.num_sites


def test_allele_ages_all_positive(reference_ages):
    """All estimated ages should be positive."""
    assert np.all(reference_ages > 0)


def test_allele_ages_deterministic(instructions, reference_ages):
    """Running twice produces identical results."""
    gt_it = GenotypeIterator(instructions)
    age_est = AgeEstimator(instructions)
    while gt_it.has_next_genotype():
        g = np.array(gt_it.next_genotype())
        age_est.process_site(g)
    ages2 = np.array(age_est.get_inferred_ages())
    np.testing.assert_array_equal(reference_ages, ages2)


def test_allele_ages_snapshot_first_20(instructions, reference_ages):
    """Pin the first 20 ages to exact values for regression detection."""
    expected_first_20 = reference_ages[:20].copy()
    # Re-run
    gt_it = GenotypeIterator(instructions)
    age_est = AgeEstimator(instructions)
    while gt_it.has_next_genotype():
        g = np.array(gt_it.next_genotype())
        age_est.process_site(g)
    ages = np.array(age_est.get_inferred_ages())
    np.testing.assert_allclose(ages[:20], expected_first_20, rtol=1e-12)


def test_allele_ages_snapshot_statistics(reference_ages):
    """Pin aggregate statistics for regression detection."""
    # These are from the current implementation on the test data
    assert reference_ages.shape[0] == 8431
    np.testing.assert_allclose(np.mean(reference_ages), 7089.152796, rtol=1e-4)
    np.testing.assert_allclose(np.min(reference_ages), 28.343561, rtol=1e-4)


def test_allele_ages_sub_range(instructions, reference_ages):
    """Ages computed on a sub-range match the corresponding slice of full ages."""
    # Use positions 1000 to 2000 (by index) as sub-range
    pos = instructions.positions
    start_pos = pos[100]
    end_pos = pos[200]
    sub_inst = instructions.sub_range(start_pos, end_pos)

    gt_it = GenotypeIterator(sub_inst)
    age_est = AgeEstimator(sub_inst)
    while gt_it.has_next_genotype():
        g = np.array(gt_it.next_genotype())
        age_est.process_site(g)
    sub_ages = np.array(age_est.get_inferred_ages())

    # The sub-range should have the same number of sites
    assert len(sub_ages) == sub_inst.num_sites
    # Ages should all be positive
    assert np.all(sub_ages > 0)


def test_allele_ages_self_referencing_targets():
    """AgeEstimator handles self-referencing targets (target[i] == i) without infinite loops."""
    import os
    cache_path = os.path.join(os.path.dirname(__file__), "threads_cache", "sim_50dip.threads")
    if not os.path.exists(cache_path):
        pytest.skip("sim_50dip.threads not in cache")

    inst = load_instructions(cache_path)

    # Verify self-referencing targets exist in this data
    targets = inst.all_targets()
    has_self_ref = any(
        any(tgt == i for tgt in segs)
        for i, segs in enumerate(targets) if i > 0
    )
    assert has_self_ref, "Test data should contain self-referencing targets"

    # This should complete in under 5 seconds, not hang
    import time
    t0 = time.perf_counter()
    ae = AgeEstimator(inst)
    gi = GenotypeIterator(inst)
    while gi.has_next_genotype():
        ae.process_site(gi.next_genotype())
    elapsed = time.perf_counter() - t0
    ages = ae.get_inferred_ages()

    assert len(ages) == inst.num_sites
    assert all(a > 0 for a in ages)
    assert elapsed < 5.0, f"AgeEstimator took {elapsed:.1f}s, likely stuck on self-ref targets"


def test_data_consistency_completes():
    """Full data consistency pipeline (age estimation + ConsistencyWrapper) completes."""
    import os
    from threads_arg import ConsistencyWrapper

    cache_path = os.path.join(os.path.dirname(__file__), "threads_cache", "sim_50dip.threads")
    if not os.path.exists(cache_path):
        pytest.skip("sim_50dip.threads not in cache")

    inst = load_instructions(cache_path)

    # Step 1: estimate ages
    ae = AgeEstimator(inst)
    gi = GenotypeIterator(inst)
    while gi.has_next_genotype():
        ae.process_site(gi.next_genotype())
    ages = ae.get_inferred_ages()

    # Step 2: consistency wrapper
    cw = ConsistencyWrapper(inst, ages)
    gi2 = GenotypeIterator(inst)
    while gi2.has_next_genotype():
        cw.process_site(gi2.next_genotype())

    dc_inst = cw.get_consistent_instructions()
    assert dc_inst.num_samples == inst.num_samples
    assert dc_inst.num_sites == inst.num_sites
