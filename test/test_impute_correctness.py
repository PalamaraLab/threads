"""
Property-based unit tests for the imputation pipeline.

Tests correctness properties of individual components rather than
relying solely on snapshot regression.
"""
import io
import numpy as np
import pytest
import types

from scipy.sparse import csr_array

from threads_arg.fwbw import (
    forwards_ls_hap,
    backwards_ls_hap,
    set_emission_probabilities,
    fwbw,
    MISSING,
)
from threads_arg.impute import (
    _active_site_arg_delta,
    MutationMap,
    MutationContainer,
    CachedPosteriorSnps,
    WriterVCF,
    Impute,
)


# ---------------------------------------------------------------------------
# Fixtures
# ---------------------------------------------------------------------------
@pytest.fixture
def small_panel():
    """5 haplotypes, 10 biallelic sites."""
    rng = np.random.default_rng(42)
    return rng.integers(0, 2, size=(10, 5)).astype(np.int8)


@pytest.fixture
def small_query(small_panel):
    """Single query haplotype, shape (1, m)."""
    rng = np.random.default_rng(99)
    m = small_panel.shape[0]
    return rng.integers(0, 2, size=(1, m)).astype(np.int8)


@pytest.fixture
def recomb_rates(small_panel):
    m = small_panel.shape[0]
    return np.full(m, 0.01)


@pytest.fixture
def mutation_rate():
    return 0.01


# ---------------------------------------------------------------------------
# set_emission_probabilities
# ---------------------------------------------------------------------------
class TestSetEmissionProbabilities:
    def test_output_shape(self, small_panel, small_query, mutation_rate):
        e = set_emission_probabilities(small_panel, small_query, mutation_rate)
        assert e.shape == (small_panel.shape[0], 2)

    def test_emissions_in_unit_interval(self, small_panel, small_query, mutation_rate):
        e = set_emission_probabilities(small_panel, small_query, mutation_rate)
        assert np.all(e >= 0)
        assert np.all(e <= 1)

    def test_match_mismatch_sum_to_one_polymorphic(self):
        """For polymorphic sites, e[:,0] + e[:,1] == 1."""
        panel = np.array([[0, 1, 0], [1, 0, 1]], dtype=np.int8)
        query = np.array([[1, 0]], dtype=np.int8)
        e = set_emission_probabilities(panel, query, 0.05)
        np.testing.assert_allclose(e[:, 0] + e[:, 1], 1.0)

    def test_invariant_site_emissions(self):
        """Invariant sites: mismatch=0, match=1 regardless of mutation rate."""
        panel = np.array([[0, 0, 0], [1, 0, 1]], dtype=np.int8)
        query = np.array([[0, 0]], dtype=np.int8)  # site 0 all-zero => invariant
        e = set_emission_probabilities(panel, query, 0.05)
        assert e[0, 0] == 0.0
        assert e[0, 1] == 1.0
        # site 1 polymorphic
        assert e[1, 0] == pytest.approx(0.05)
        assert e[1, 1] == pytest.approx(0.95)

    def test_mutation_rate_none(self, small_panel, small_query):
        """Auto-computed mutation rate should still produce valid emissions."""
        e = set_emission_probabilities(small_panel, small_query, None)
        assert np.all(e >= 0)
        assert np.all(e <= 1)


# ---------------------------------------------------------------------------
# Forward-backward algorithm
# ---------------------------------------------------------------------------
class TestForwardBackward:
    def test_forward_scaling_factors_positive(self, small_panel, small_query,
                                               recomb_rates, mutation_rate):
        e = set_emission_probabilities(small_panel, small_query, mutation_rate)
        F, c = forwards_ls_hap(small_panel.astype(np.float64), small_query.ravel().astype(np.float64), e, recomb_rates)
        assert np.all(c > 0)

    def test_forward_values_non_negative(self, small_panel, small_query,
                                          recomb_rates, mutation_rate):
        e = set_emission_probabilities(small_panel, small_query, mutation_rate)
        F, c = forwards_ls_hap(small_panel.astype(np.float64), small_query.ravel().astype(np.float64), e, recomb_rates)
        assert np.all(F >= 0)

    def test_forward_normalized_rows_sum_to_one(self, small_panel, small_query,
                                                  recomb_rates, mutation_rate):
        e = set_emission_probabilities(small_panel, small_query, mutation_rate)
        F, c = forwards_ls_hap(small_panel.astype(np.float64), small_query.ravel().astype(np.float64), e, recomb_rates)
        np.testing.assert_allclose(F.sum(axis=1), 1.0, atol=1e-10)

    def test_backward_last_row_all_ones(self, small_panel, small_query,
                                         recomb_rates, mutation_rate):
        m, n = small_panel.shape
        e = set_emission_probabilities(small_panel, small_query, mutation_rate)
        F, c = forwards_ls_hap(small_panel.astype(np.float64), small_query.ravel().astype(np.float64), e, recomb_rates)
        B = backwards_ls_hap(small_panel.astype(np.float64), small_query.ravel().astype(np.float64), e, c, recomb_rates)
        np.testing.assert_array_equal(B[-1], np.ones(n))

    def test_backward_values_non_negative(self, small_panel, small_query,
                                           recomb_rates, mutation_rate):
        e = set_emission_probabilities(small_panel, small_query, mutation_rate)
        F, c = forwards_ls_hap(small_panel.astype(np.float64), small_query.ravel().astype(np.float64), e, recomb_rates)
        B = backwards_ls_hap(small_panel.astype(np.float64), small_query.ravel().astype(np.float64), e, c, recomb_rates)
        assert np.all(B >= 0)

    def test_posterior_shape(self, small_panel, small_query, recomb_rates, mutation_rate):
        posterior = fwbw(small_panel, small_query, recomb_rates, mutation_rate)
        assert posterior.shape == small_panel.shape

    def test_posterior_non_negative(self, small_panel, small_query,
                                     recomb_rates, mutation_rate):
        posterior = fwbw(small_panel, small_query, recomb_rates, mutation_rate)
        assert np.all(posterior >= -1e-15)

    def test_posterior_rows_sum_to_one(self, small_panel, small_query,
                                        recomb_rates, mutation_rate):
        posterior = fwbw(small_panel, small_query, recomb_rates, mutation_rate)
        np.testing.assert_allclose(posterior.sum(axis=1), 1.0, atol=1e-10)

    def test_identical_query_concentrates_posterior(self):
        """Posterior should concentrate on haplotypes identical to query."""
        panel = np.array([[0, 1, 0, 1],
                          [1, 0, 1, 0],
                          [0, 1, 0, 1],
                          [1, 0, 1, 0]], dtype=np.int8)
        # h0 == h2 == [0,1,0,1]; h1 == h3 == [1,0,1,0]
        query = panel[:, 0:1].T  # (1, m): matches h0 and h2
        recomb = np.full(4, 0.001)
        posterior = fwbw(panel, query, recomb, 0.001)
        match_weight = posterior[:, [0, 2]].sum(axis=1)
        other_weight = posterior[:, [1, 3]].sum(axis=1)
        assert np.all(match_weight > other_weight)

    def test_missing_data_still_valid(self):
        """MISSING sites should not break posterior properties."""
        panel = np.array([[0, 1, 0], [1, 0, 1], [0, 1, 0]], dtype=np.int8)
        query = np.array([[MISSING, 1, 0]], dtype=np.int8)
        recomb = np.full(3, 0.01)
        posterior = fwbw(panel, query, recomb, 0.01)
        assert np.all(posterior >= -1e-15)
        np.testing.assert_allclose(posterior.sum(axis=1), 1.0, atol=1e-10)

    def test_different_panel_sizes(self):
        """Properties hold for various n and m."""
        for m, n in [(3, 2), (5, 10), (20, 3), (50, 50)]:
            rng = np.random.default_rng(m * 100 + n)
            panel = rng.integers(0, 2, size=(m, n)).astype(np.int8)
            query = rng.integers(0, 2, size=(1, m)).astype(np.int8)
            recomb = np.full(m, 0.01)
            posterior = fwbw(panel, query, recomb, 0.01)
            assert posterior.shape == (m, n)
            assert np.all(posterior >= -1e-15)
            np.testing.assert_allclose(posterior.sum(axis=1), 1.0, atol=1e-10)


# ---------------------------------------------------------------------------
# _sparsify_posterior
# ---------------------------------------------------------------------------
class TestSparsifyPosterior:
    def _mock(self, num_snps, num_samples_panel):
        return types.SimpleNamespace(
            num_snps=num_snps,
            num_samples_panel=num_samples_panel,
        )

    def test_output_shape(self):
        mock = self._mock(10, 100)
        posterior = np.random.default_rng(0).random((10, 5)) + 0.01
        matched = np.array([3, 17, 42, 55, 88])
        result = Impute._sparsify_posterior(mock, posterior, matched)
        assert result.shape == (10, 100)

    def test_rows_sum_to_one(self):
        mock = self._mock(10, 100)
        posterior = np.random.default_rng(0).random((10, 5)) + 0.01
        matched = np.array([3, 17, 42, 55, 88])
        result = Impute._sparsify_posterior(mock, posterior, matched)
        row_sums = np.asarray(result.sum(axis=1)).ravel()
        np.testing.assert_allclose(row_sums, 1.0, atol=1e-10)

    def test_no_negative_values(self):
        mock = self._mock(10, 50)
        posterior = np.random.default_rng(0).random((10, 5)) + 0.01
        matched = np.arange(5)
        result = Impute._sparsify_posterior(mock, posterior, matched)
        assert np.all(result.data >= 0)

    def test_columns_in_matched_set(self):
        mock = self._mock(5, 100)
        posterior = np.ones((5, 3)) * 0.5
        matched = np.array([10, 50, 90])
        result = Impute._sparsify_posterior(mock, posterior, matched)
        nz_cols = set(result.nonzero()[1])
        assert nz_cols.issubset({10, 50, 90})

    def test_small_values_thresholded(self):
        """Entries <= 1/num_samples_panel should be dropped."""
        mock = self._mock(3, 10)
        # Row 0: one value well above threshold, rest tiny
        posterior = np.array([[1.0, 0.001, 0.001],
                              [0.5, 0.4, 0.1],
                              [0.9, 0.05, 0.05]])
        matched = np.array([0, 1, 2])
        result = Impute._sparsify_posterior(mock, posterior, matched)
        # threshold = 1/10 = 0.1; row 0 col 1,2 (0.001) should be gone
        row0 = np.asarray(result[0].todense()).ravel()
        assert row0[1] == 0.0
        assert row0[2] == 0.0
        assert row0[0] > 0.0


# ---------------------------------------------------------------------------
# MutationMap
# ---------------------------------------------------------------------------
class TestMutationMap:
    def test_unmapped(self):
        mm = MutationMap("snp1", 0, "NaN")
        assert not mm.is_mapped()
        assert not mm.flipped

    def test_single_uniquely_mapped(self):
        mm = MutationMap("snp1", 0, "-1,100.0,200.0")
        assert mm.is_mapped()
        assert mm.uniquely_mapped
        assert mm.get_boundaries(-1) == (100.0, 200.0)
        # Uniquely mapped: any query returns the -1 boundaries
        assert mm.get_boundaries(42) == (100.0, 200.0)

    def test_multiple_mappings(self):
        mm = MutationMap("snp1", 0, "5,100.0,200.0;10,300.0,400.0")
        assert mm.is_mapped()
        assert not mm.uniquely_mapped
        assert mm.get_boundaries(5) == (100.0, 200.0)
        assert mm.get_boundaries(10) == (300.0, 400.0)

    def test_dotted_sample_ids(self):
        """'5.10.15' expands to three separate carrier IDs."""
        mm = MutationMap("snp1", 0, "5.10.15,100.0,200.0")
        assert mm.is_carrier(5)
        assert mm.is_carrier(10)
        assert mm.is_carrier(15)
        assert not mm.is_carrier(20)

    def test_flipped_flag(self):
        mm = MutationMap("snp1", 1, "-1,100.0,200.0")
        assert mm.flipped

    def test_multi_group_carriers(self):
        mm = MutationMap("snp1", 0, "1.2,10.0,20.0;3,30.0,40.0")
        assert mm.is_carrier(1)
        assert mm.is_carrier(2)
        assert mm.is_carrier(3)
        assert mm.get_boundaries(1) == (10.0, 20.0)
        assert mm.get_boundaries(3) == (30.0, 40.0)


# ---------------------------------------------------------------------------
# MutationContainer
# ---------------------------------------------------------------------------
class TestMutationContainer:
    def test_load_and_lookup(self, tmp_path):
        mut_file = tmp_path / "test.mut"
        mut_file.write_text(
            "'snp1'\t100\t0\t-1,50.0,150.0\n"
            "'snp2'\t200\t0\tNaN\n"
            "'snp3'\t300\t1\t5.10,100.0,200.0;15,300.0,400.0\n"
        )
        mc = MutationContainer(str(mut_file))
        assert mc.is_mapped("'snp1'")
        assert not mc.is_mapped("'snp2'")
        assert mc.is_mapped("'snp3'")
        assert not mc.is_mapped("nonexistent")

    def test_get_mapping(self, tmp_path):
        mut_file = tmp_path / "test.mut"
        mut_file.write_text("'snp1'\t100\t0\t-1,50.0,150.0\n")
        mc = MutationContainer(str(mut_file))
        mm = mc.get_mapping("'snp1'")
        assert mm.uniquely_mapped
        assert mm.get_boundaries(-1) == (50.0, 150.0)


# ---------------------------------------------------------------------------
# _active_site_arg_delta
# ---------------------------------------------------------------------------
class TestActiveSiteArgDelta:
    def _seg(self, seg_start, ids, ages):
        return types.SimpleNamespace(seg_start=seg_start, ids=ids, ages=ages)

    def _rec(self, pos):
        return types.SimpleNamespace(pos=pos)

    def test_zero_delta_no_carriers_in_active(self):
        posterior = np.array([0.5, 0.3, 0.2])
        active_indexes = {0: 0, 1: 1, 2: 2}
        thread = [self._seg(0, [0, 1, 2], [100.0, 200.0, 300.0])]
        mm = MutationMap("snp", 0, "99,50.0,150.0")  # carrier 99 not active
        carriers = {99}
        delta = _active_site_arg_delta(
            posterior, active_indexes, thread, mm, carriers, self._rec(50))
        assert delta == 0.0

    def test_zero_delta_zero_posterior(self):
        posterior = np.array([0.0, 0.5, 0.5])
        active_indexes = {0: 0, 1: 1, 2: 2}
        thread = [self._seg(0, [0, 1, 2], [100.0, 200.0, 300.0])]
        mm = MutationMap("snp", 0, "0,50.0,150.0")  # carrier 0, posterior=0
        carriers = {0}
        delta = _active_site_arg_delta(
            posterior, active_indexes, thread, mm, carriers, self._rec(50))
        assert delta == 0.0

    def test_delta_sign_is_negative(self):
        """arg_prob < 1 => (arg_prob - 1) < 0 => delta < 0 when posterior > 0."""
        posterior = np.array([0.8, 0.1, 0.1])
        active_indexes = {0: 0}
        thread = [self._seg(0, [0], [1000.0])]
        # Use -1 sentinel for uniquely mapped (matches real .mut format)
        mm = MutationMap("snp", 0, "-1,50.0,150.0")
        carriers = {0}
        delta = _active_site_arg_delta(
            posterior, active_indexes, thread, mm, carriers, self._rec(50))
        assert delta < 0

    def test_erlang2_probability_in_unit_interval(self):
        """The Erlang(2) CDF must be in [0, 1] for all positive inputs."""
        for height in [10.0, 100.0, 1000.0, 10000.0]:
            for mut_h in [1.0, 50.0, 500.0]:
                lam = 2.0 / height
                lam_mut = lam * mut_h
                arg_prob = 1 - np.exp(-lam_mut) * (1 + lam_mut)
                assert 0 <= arg_prob <= 1, f"height={height}, mut_h={mut_h}"

    def test_segment_selection_by_position(self):
        """Position determines which segment's IDs are used."""
        posterior = np.array([0.5, 0.5])
        active_indexes = {0: 0, 1: 1}
        seg1 = self._seg(0, [0], [100.0])
        seg2 = self._seg(500, [1], [100.0])
        thread = [seg1, seg2]
        # Non-uniquely mapped: both sample 0 and 1 have explicit boundaries
        mm = MutationMap("snp", 0, "0,50.0,150.0;1,50.0,150.0")
        carriers = {0, 1}

        # pos=250 => seg1 (carrier 0)
        d1 = _active_site_arg_delta(
            posterior, active_indexes, thread, mm, carriers, self._rec(250))
        # pos=750 => seg2 (carrier 1)
        d2 = _active_site_arg_delta(
            posterior, active_indexes, thread, mm, carriers, self._rec(750))
        assert d1 != 0.0
        assert d2 != 0.0

    def test_young_mutation_stronger_correction(self):
        """Young mutation (low mut_height) relative to coalescence => larger |delta|.

        When mut_height << coalescence height, the Erlang(2) prob is small,
        so (arg_prob - 1) ≈ -1 and the correction removes most posterior weight.
        When mut_height ≈ height, arg_prob ≈ 1 and little correction occurs.
        """
        posterior = np.array([0.5])
        active_indexes = {0: 0}
        thread = [self._seg(0, [0], [500.0])]
        carriers = {0}

        mm_young = MutationMap("snp", 0, "-1,10.0,20.0")    # mut_height ~15
        mm_old = MutationMap("snp", 0, "-1,400.0,500.0")    # mut_height ~450

        d_young = _active_site_arg_delta(
            posterior, active_indexes, thread, mm_young, carriers, self._rec(50))
        d_old = _active_site_arg_delta(
            posterior, active_indexes, thread, mm_old, carriers, self._rec(50))
        assert abs(d_young) > abs(d_old)


# ---------------------------------------------------------------------------
# CachedPosteriorSnps
# ---------------------------------------------------------------------------
class TestCachedPosteriorSnps:
    def _make_posteriors(self, n_targets=3, n_snps=5, n_panel=20, seed=42):
        rng = np.random.default_rng(seed)
        posteriors = []
        for _ in range(n_targets):
            dense = rng.random((n_snps, n_panel))
            dense[dense < 0.7] = 0
            # Guarantee at least one nonzero per row
            for i in range(n_snps):
                if dense[i].sum() == 0:
                    dense[i, rng.integers(n_panel)] = 1.0
            row_sums = dense.sum(axis=1, keepdims=True)
            dense /= row_sums
            posteriors.append(csr_array(dense))
        return posteriors

    def test_output_shape(self):
        posteriors = self._make_posteriors(n_targets=3, n_panel=20)
        cache = CachedPosteriorSnps(posteriors)
        assert cache[0].shape == (3, 20)

    def test_rows_normalized(self):
        posteriors = self._make_posteriors()
        cache = CachedPosteriorSnps(posteriors)
        result = cache[0]
        np.testing.assert_allclose(result.sum(axis=1), 1.0, atol=1e-10)

    def test_cache_hit_same_object(self):
        posteriors = self._make_posteriors()
        cache = CachedPosteriorSnps(posteriors)
        r1 = cache[0]
        r2 = cache[0]
        assert r1 is r2

    def test_cache_eviction(self):
        posteriors = self._make_posteriors()
        cache = CachedPosteriorSnps(posteriors, max_size=2)
        _ = cache[0]
        _ = cache[1]
        assert 0 in cache.posteriors_by_snp_idx
        _ = cache[2]  # evicts 0
        assert 0 not in cache.posteriors_by_snp_idx
        assert 1 in cache.posteriors_by_snp_idx
        assert 2 in cache.posteriors_by_snp_idx

    def test_negative_index(self):
        posteriors = self._make_posteriors(n_snps=5)
        cache = CachedPosteriorSnps(posteriors)
        r_neg = cache[-1]
        cache2 = CachedPosteriorSnps(posteriors)
        r_pos = cache2[4]
        np.testing.assert_array_equal(r_neg, r_pos)


# ---------------------------------------------------------------------------
# WriterVCF format correctness
# ---------------------------------------------------------------------------
class TestWriterVCF:
    def _rec(self, pos=100, ref="A", alt=["T"], af=0.1, id="snp1"):
        return types.SimpleNamespace(pos=pos, ref=ref, alt=alt, af=af, id=id)

    def _write_line(self, genotypes, imputed=True):
        buf = io.StringIO()
        writer = WriterVCF(None)
        writer.file = buf
        writer.write_site(genotypes, self._rec(), imputed, "1")
        return buf.getvalue().strip()

    def test_gt_field_is_binary(self):
        line = self._write_line(np.array([0.0, 1.0, 0.3, 0.7]))
        for sample in line.split("\t")[9:]:
            gt = sample.split(":")[0]
            h1, h2 = gt.split("|")
            assert h1 in ("0", "1")
            assert h2 in ("0", "1")

    def test_dosage_equals_hap_sum(self):
        genotypes = np.array([0.2, 0.8, 0.0, 1.0])
        line = self._write_line(genotypes)
        for i, sample in enumerate(line.split("\t")[9:]):
            ds = float(sample.split(":")[1])
            expected = round(genotypes[2*i] + genotypes[2*i+1], 3)
            assert abs(ds - expected) < 1e-6

    def test_imp_flag_present_when_imputed(self):
        assert "IMP;" in self._write_line(np.array([0.0, 1.0]), imputed=True)

    def test_imp_flag_absent_when_not_imputed(self):
        assert "IMP;" not in self._write_line(np.array([0.0, 1.0]), imputed=False)

    def test_trailing_zero_stripping(self):
        # hap1=0.0, hap2=0.3 => dosage=0.3
        line = self._write_line(np.array([0.0, 0.3]))
        ds_str = line.split("\t")[9].split(":")[1]
        assert ds_str == "0.3"

    def test_integer_dosage_no_decimal(self):
        # hap1=0.0, hap2=1.0 => dosage=1.0 => "1"
        line = self._write_line(np.array([0.0, 1.0]))
        ds_str = line.split("\t")[9].split(":")[1]
        assert ds_str == "1"

    def test_zero_dosage_format(self):
        line = self._write_line(np.array([0.0, 0.0]))
        ds_str = line.split("\t")[9].split(":")[1]
        assert ds_str == "0"


# ---------------------------------------------------------------------------
# Dosage bounds (integration property)
# ---------------------------------------------------------------------------
class TestDosageBounds:
    def test_haploid_dosages_in_unit_interval(self):
        """Posterior rows sum to 1, so any column-subset sum is in [0, 1]."""
        rng = np.random.default_rng(42)
        for _ in range(5):
            m = rng.integers(5, 30)
            n = rng.integers(3, 20)
            panel = rng.integers(0, 2, size=(m, n)).astype(np.int8)
            query = rng.integers(0, 2, size=(1, m)).astype(np.int8)
            recomb = np.full(m, 0.01)
            posterior = fwbw(panel, query, recomb, 0.01)
            # Any boolean mask over columns sums to <= 1 per row
            mask = rng.choice([True, False], size=n)
            if not mask.any():
                mask[0] = True
            subset_sums = posterior[:, mask].sum(axis=1)
            assert np.all(subset_sums >= -1e-10)
            assert np.all(subset_sums <= 1 + 1e-10)
