"""
Tests for the `threads infer` command: utility functions, Matcher, ThreadsLowMem,
ThreadingInstructions, ConsistencyWrapper, and serialization round-trips.
"""
import pickle
import tempfile

import numpy as np
import pytest

from pathlib import Path

TEST_DATA = Path(__file__).parent / "data"
PANEL_PGEN = str(TEST_DATA / "panel.pgen")
PANEL_PVAR = str(TEST_DATA / "panel.pvar")
PANEL_PSAM = str(TEST_DATA / "panel.psam")
GMAP = str(TEST_DATA / "gmap_02.map")
DEMO = str(TEST_DATA / "CEU_unscaled.demo")
THREADS_FIT = str(TEST_DATA / "expected_infer_fit_to_data_snapshot.threads")
THREADS_NO_FIT = str(TEST_DATA / "expected_infer_snapshot.threads")


# ===================================================================
# Utility functions: split_list
# ===================================================================
class TestSplitList:
    def test_even_split(self):
        from threads_arg.utils import split_list
        result = split_list([1, 2, 3, 4, 5, 6], 3)
        assert len(result) == 3
        assert sum(len(s) for s in result) == 6
        # Flattened should be original
        assert [x for s in result for x in s] == [1, 2, 3, 4, 5, 6]

    def test_uneven_split(self):
        from threads_arg.utils import split_list
        result = split_list([1, 2, 3, 4, 5], 3)
        assert len(result) == 3
        assert sum(len(s) for s in result) == 5
        # First chunks get the extra elements
        sizes = [len(s) for s in result]
        assert sizes == [2, 2, 1]

    def test_n_equals_length(self):
        from threads_arg.utils import split_list
        result = split_list([1, 2, 3], 3)
        assert result == [[1], [2], [3]]

    def test_n_greater_than_length(self):
        from threads_arg.utils import split_list
        result = split_list([1, 2], 5)
        assert len(result) == 5
        non_empty = [s for s in result if len(s) > 0]
        assert [x for s in non_empty for x in s] == [1, 2]

    def test_single_chunk(self):
        from threads_arg.utils import split_list
        result = split_list([1, 2, 3, 4], 1)
        assert result == [[1, 2, 3, 4]]

    def test_empty_list(self):
        from threads_arg.utils import split_list
        result = split_list([], 3)
        assert len(result) == 3
        assert all(len(s) == 0 for s in result)


# ===================================================================
# Utility functions: parse_demography
# ===================================================================
class TestParseDemography:
    def test_loads_ceu_demo(self):
        from threads_arg.utils import parse_demography
        times, sizes = parse_demography(DEMO)
        assert len(times) == len(sizes)
        assert len(times) > 0
        assert times[0] == 0.0
        assert all(t >= 0 for t in times)
        assert all(s > 0 for s in sizes)

    def test_times_monotonic(self):
        from threads_arg.utils import parse_demography
        times, _ = parse_demography(DEMO)
        for i in range(1, len(times)):
            assert times[i] > times[i - 1]

    def test_custom_demo_file(self):
        from threads_arg.utils import parse_demography
        with tempfile.NamedTemporaryFile(mode="w", suffix=".demo", delete=False) as f:
            f.write("0.0\t10000\n100.0\t20000\n500.0\t5000\n")
            f.flush()
            times, sizes = parse_demography(f.name)
        assert times == [0.0, 100.0, 500.0]
        assert sizes == [10000.0, 20000.0, 5000.0]


# ===================================================================
# Utility functions: recombination maps
# ===================================================================
class TestRecombination:
    def test_constant_recombination_shape(self):
        from threads_arg.utils import make_constant_recombination_from_pgen
        cm, phys = make_constant_recombination_from_pgen(PANEL_PGEN, 1.3e-8)
        assert len(cm) == len(phys)
        assert len(cm) > 0

    def test_constant_recombination_monotonic(self):
        from threads_arg.utils import make_constant_recombination_from_pgen
        cm, phys = make_constant_recombination_from_pgen(PANEL_PGEN, 1.3e-8)
        assert np.all(np.diff(cm) > 0), "cM positions must be strictly increasing"
        assert np.all(np.diff(phys) > 0), "physical positions must be strictly increasing"

    def test_map_recombination_shape(self):
        from threads_arg.utils import make_recombination_from_map_and_pgen
        # gmap_02.map uses chr 20, pass None to skip chromosome check
        cm, phys = make_recombination_from_map_and_pgen(GMAP, PANEL_PGEN, None)
        assert len(cm) == len(phys)
        assert len(cm) > 0

    def test_map_recombination_monotonic(self):
        from threads_arg.utils import make_recombination_from_map_and_pgen
        cm, phys = make_recombination_from_map_and_pgen(GMAP, PANEL_PGEN, None)
        assert np.all(np.diff(cm) > 0)


# ===================================================================
# Utility functions: read_all_genotypes, read_sample_names, etc.
# ===================================================================
class TestPgenReading:
    def test_read_all_genotypes_shape(self):
        from threads_arg.utils import read_all_genotypes
        gt = read_all_genotypes(PANEL_PGEN)
        assert gt.ndim == 2
        assert gt.dtype == np.int32
        # 500 diploid samples = 1000 haploids
        assert gt.shape[1] == 1000

    def test_read_all_genotypes_biallelic(self):
        from threads_arg.utils import read_all_genotypes
        gt = read_all_genotypes(PANEL_PGEN)
        assert set(np.unique(gt)).issubset({0, 1})

    def test_read_sample_names(self):
        from threads_arg.utils import read_sample_names
        names = read_sample_names(PANEL_PGEN)
        assert len(names) == 500
        assert names[0] == "tsk_0"
        assert names[-1] == "tsk_499"

    def test_read_positions_and_ids(self):
        from threads_arg.utils import read_positions_and_ids
        positions, ids = read_positions_and_ids(PANEL_PGEN)
        assert len(positions) == len(ids)
        assert len(positions) > 0
        assert all(isinstance(p, int) for p in positions)
        # Positions should be sorted
        assert positions == sorted(positions)

    def test_read_variant_metadata_columns(self):
        from threads_arg.utils import read_variant_metadata
        df = read_variant_metadata(PANEL_PGEN)
        for col in ["CHROM", "POS", "ID", "REF", "ALT"]:
            assert col in df.columns
        assert len(df) > 0


# ===================================================================
# Matcher: PBWT haplotype matching
# ===================================================================
class TestMatcher:
    @pytest.fixture(scope="class")
    def matcher_setup(self):
        from threads_arg import Matcher
        from threads_arg.utils import (
            make_recombination_from_map_and_pgen,
            read_all_genotypes,
        )
        cm, phys = make_recombination_from_map_and_pgen(GMAP, PANEL_PGEN, None)
        gt = read_all_genotypes(PANEL_PGEN)
        n_haps = gt.shape[1]

        # Filter singletons
        ac = gt.sum(axis=1)
        ac_mask = (ac > 1) & (ac < n_haps)

        matcher = Matcher(n_haps, cm[ac_mask], 0.01, 0.5, 4, 4)
        if hasattr(matcher, "process_all_sites_numpy"):
            matcher.process_all_sites_numpy(gt[ac_mask])
        else:
            for g in gt[ac_mask]:
                matcher.process_site(g)
        matcher.propagate_adjacent_matches()
        return matcher, n_haps

    def test_num_samples(self, matcher_setup):
        matcher, n_haps = matcher_setup
        assert matcher.num_samples == n_haps

    def test_num_sites(self, matcher_setup):
        matcher, _ = matcher_setup
        assert matcher.num_sites > 0

    def test_get_matches_returns_list(self, matcher_setup):
        matcher, _ = matcher_setup
        # get_matches() takes no args (returns all match groups)
        matches = matcher.get_matches()
        assert isinstance(matches, list)
        assert len(matches) > 0

    def test_match_groups_have_candidates(self, matcher_setup):
        matcher, _ = matcher_setup
        matches = matcher.get_matches()
        for mg in matches[:10]:  # spot-check first 10
            assert len(mg.match_candidates) > 0

    def test_cm_positions_length(self, matcher_setup):
        matcher, _ = matcher_setup
        cm_pos = matcher.cm_positions()
        assert len(cm_pos) > 0

    def test_serializable_matches_shape(self, matcher_setup):
        matcher, _ = matcher_setup
        sample_ids = [0, 1, 2]
        s_matches = matcher.serializable_matches(sample_ids)
        assert len(s_matches) > 0


# ===================================================================
# ThreadsLowMem: HMM inference engine
# ===================================================================
class TestThreadsLowMem:
    @pytest.fixture(scope="class")
    def small_inference(self):
        """Run a small inference on 10 haploids to test the C++ engine."""
        from threads_arg import ThreadsLowMem, Matcher
        from threads_arg.utils import (
            make_recombination_from_map_and_pgen,
            read_all_genotypes,
            parse_demography,
        )
        cm, phys = make_recombination_from_map_and_pgen(GMAP, PANEL_PGEN, None)
        gt = read_all_genotypes(PANEL_PGEN)
        n_haps = gt.shape[1]
        ne_times, ne = parse_demography(DEMO)

        ac = gt.sum(axis=1)
        ac_mask = (ac > 1) & (ac < n_haps)

        matcher = Matcher(n_haps, cm[ac_mask], 0.01, 0.5, 4, 4)
        if hasattr(matcher, "process_all_sites_numpy"):
            matcher.process_all_sites_numpy(gt[ac_mask])
        else:
            for g in gt[ac_mask]:
                matcher.process_site(g)
        matcher.propagate_adjacent_matches()

        # Infer just 10 haploids
        target_ids = list(range(10))
        s_matches = matcher.serializable_matches(target_ids)
        cm_pos = matcher.cm_positions()

        tlm = ThreadsLowMem(target_ids, phys, cm, ne, ne_times, 1.4e-8, False)
        tlm.initialize_viterbi(s_matches, cm_pos)

        if hasattr(tlm, "process_all_sites_viterbi_numpy"):
            tlm.process_all_sites_viterbi_numpy(gt)
        else:
            for g in gt:
                tlm.process_site_viterbi(g)

        tlm.prune()
        tlm.traceback()

        if hasattr(tlm, "process_all_sites_hets_numpy"):
            tlm.process_all_sites_hets_numpy(gt)
        else:
            for g in gt:
                tlm.process_site_hets(g)

        tlm.date_segments()
        return tlm, target_ids, phys

    def test_serialize_paths_length(self, small_inference):
        tlm, target_ids, _ = small_inference
        seg_starts, match_ids, heights, hetsites = tlm.serialize_paths()
        assert len(seg_starts) == len(target_ids)
        assert len(match_ids) == len(target_ids)
        assert len(heights) == len(target_ids)
        assert len(hetsites) == len(target_ids)

    def test_segment_starts_mostly_non_empty(self, small_inference):
        """Most samples should have at least one segment (sample 0 may be empty as it's the reference)."""
        tlm, target_ids, _ = small_inference
        seg_starts, _, _, _ = tlm.serialize_paths()
        non_empty = sum(1 for ss in seg_starts if len(ss) > 0)
        # At least all non-reference samples should have segments
        assert non_empty >= len(target_ids) - 1

    def test_heights_positive(self, small_inference):
        tlm, _, _ = small_inference
        _, _, heights, _ = tlm.serialize_paths()
        for h_list in heights:
            for h in h_list:
                assert h >= 0, f"negative height {h}"

    def test_match_ids_valid(self, small_inference):
        tlm, target_ids, _ = small_inference
        _, match_ids, _, _ = tlm.serialize_paths()
        n_haps = 1000  # panel.pgen
        for mi_list in match_ids:
            for mi in mi_list:
                assert 0 <= mi < n_haps, f"match_id {mi} out of range"

    def test_hetsites_within_bounds(self, small_inference):
        tlm, _, phys = small_inference
        _, _, _, hetsites = tlm.serialize_paths()
        n_sites = len(phys)
        for hs_list in hetsites:
            for hs in hs_list:
                assert 0 <= hs < n_sites, f"hetsite index {hs} out of range"


# ===================================================================
# ThreadingInstructions: construction, accessors, sub_range
# ===================================================================
class TestThreadingInstructions:
    def test_construction_from_lists(self):
        from threads_arg import ThreadingInstructions
        starts = [[0, 50], [0, 30, 70]]
        tmrcas = [[100.0, 200.0], [150.0, 250.0, 300.0]]
        targets = [[1, 0], [0, 1, 0]]
        mismatches = [[2, 4], [1, 3, 5]]
        positions = [10, 20, 30, 40, 50, 60, 70, 80]
        ti = ThreadingInstructions(starts, tmrcas, targets, mismatches, positions, 0, 100)

        assert ti.num_samples == 2
        assert ti.num_sites == 8
        assert ti.start == 0
        assert ti.end == 100
        assert ti.positions == positions

    def test_all_starts_accessor(self):
        from threads_arg import ThreadingInstructions
        starts = [[0, 50], [0, 30]]
        ti = ThreadingInstructions(starts, [[1.0, 2.0], [3.0, 4.0]],
                                   [[1, 0], [0, 1]], [[], []], [10, 20, 30], 0, 50)
        assert ti.all_starts() == starts

    def test_all_tmrcas_accessor(self):
        from threads_arg import ThreadingInstructions
        tmrcas = [[100.0, 200.0], [300.0, 400.0]]
        ti = ThreadingInstructions([[0, 50], [0, 30]], tmrcas,
                                   [[1, 0], [0, 1]], [[], []], [10, 20, 30], 0, 50)
        assert ti.all_tmrcas() == tmrcas

    def test_sub_range(self):
        from threads_arg import ThreadingInstructions
        positions = [10, 20, 30, 40, 50]
        ti = ThreadingInstructions(
            [[0], [0]], [[100.0], [200.0]], [[1], [0]],
            [[1, 3], [0, 2, 4]], positions, 0, 60
        )
        sub = ti.sub_range(15, 45)
        # sub_range should reduce the number of sites
        assert sub.num_sites <= ti.num_sites
        assert sub.num_samples == ti.num_samples
        # start/end are snapped to actual variant positions within range
        assert sub.start >= 15
        assert sub.end <= 60  # original end

    def test_pickle_roundtrip(self):
        from threads_arg import ThreadingInstructions
        ti = ThreadingInstructions(
            [[0, 50]], [[100.0, 200.0]], [[1, 0]], [[2, 4]],
            [10, 20, 30, 40, 50], 10, 50
        )
        restored = pickle.loads(pickle.dumps(ti))
        assert restored.all_starts() == ti.all_starts()
        assert restored.all_tmrcas() == ti.all_tmrcas()
        assert restored.all_targets() == ti.all_targets()
        assert restored.all_mismatches() == ti.all_mismatches()
        assert restored.positions == ti.positions
        assert restored.start == ti.start
        assert restored.end == ti.end

    def test_load_from_threads_file(self):
        from threads_arg.serialization import load_instructions
        ti = load_instructions(THREADS_NO_FIT)
        assert ti.num_samples == 1000
        assert ti.num_sites == 8431
        assert len(ti.positions) == 8431


# ===================================================================
# ViterbiPath: construction and accessors
# ===================================================================
class TestViterbiPath:
    def test_construct_empty(self):
        from threads_arg import ViterbiPath
        vp = ViterbiPath(42)
        assert vp.size() == 0

    def test_construct_with_data(self):
        from threads_arg import ViterbiPath
        vp = ViterbiPath(0, [0, 100, 500], [5, 10, 3], [100.0, 200.0, 150.0], [2, 7])
        assert vp.size() == 3
        assert vp.segment_starts == [0, 100, 500]
        assert vp.sample_ids == [5, 10, 3]
        assert vp.het_sites == [2, 7]

    def test_heights_accessible(self):
        from threads_arg import ViterbiPath
        vp = ViterbiPath(0, [0, 100], [5, 10], [123.4, 567.8], [])
        assert vp.heights == pytest.approx([123.4, 567.8])


# ===================================================================
# ConsistencyWrapper: fit-to-data post-processing
# ===================================================================
class TestConsistencyWrapper:
    def test_consistent_instructions_shape(self):
        from threads_arg import ConsistencyWrapper, GenotypeIterator
        from threads_arg.serialization import load_instructions

        ti = load_instructions(THREADS_FIT)
        # The fit_to_data .threads has embedded allele ages
        import h5py
        with h5py.File(THREADS_FIT, "r") as f:
            ages = f["allele_ages"][:].tolist()

        cw = ConsistencyWrapper(ti, ages)
        gt_it = GenotypeIterator(ti)
        while gt_it.has_next_genotype():
            g = np.array(gt_it.next_genotype())
            cw.process_site(g)

        consistent = cw.get_consistent_instructions()
        assert consistent.num_samples == ti.num_samples
        assert consistent.num_sites == ti.num_sites

    def test_consistent_preserves_genotypes(self):
        """Consistency wrapper should not change the reconstructed genotypes."""
        from threads_arg import ConsistencyWrapper, GenotypeIterator
        from threads_arg.serialization import load_instructions

        ti = load_instructions(THREADS_FIT)
        import h5py
        with h5py.File(THREADS_FIT, "r") as f:
            ages = f["allele_ages"][:].tolist()

        cw = ConsistencyWrapper(ti, ages)
        gt_it = GenotypeIterator(ti)
        original_genotypes = []
        while gt_it.has_next_genotype():
            g = np.array(gt_it.next_genotype())
            original_genotypes.append(g.copy())
            cw.process_site(g)

        consistent = cw.get_consistent_instructions()
        gt_it2 = GenotypeIterator(consistent)
        for i, orig_g in enumerate(original_genotypes):
            new_g = np.array(gt_it2.next_genotype())
            np.testing.assert_array_equal(
                new_g, orig_g,
                err_msg=f"genotype changed at site {i}"
            )


# ===================================================================
# Serialization: round-trip .threads write/read
# ===================================================================
class TestSerialization:
    def test_serialize_load_roundtrip(self):
        from threads_arg import ThreadingInstructions
        from threads_arg.serialization import serialize_instructions, load_instructions

        starts = [[0, 50], [0, 30, 70]]
        tmrcas = [[100.0, 200.0], [150.0, 250.0, 300.0]]
        targets = [[1, 0], [0, 1, 0]]
        mismatches = [[2, 4], [1, 3, 5]]
        positions = [10, 20, 30, 40, 50, 60, 70, 80]
        ti = ThreadingInstructions(starts, tmrcas, targets, mismatches, positions, 0, 100)

        with tempfile.NamedTemporaryFile(suffix=".threads", delete=False) as f:
            out_path = f.name

        serialize_instructions(ti, out_path)
        loaded = load_instructions(out_path)

        assert loaded.num_samples == ti.num_samples
        assert loaded.num_sites == ti.num_sites
        assert loaded.positions == ti.positions
        assert loaded.start == ti.start
        assert loaded.end == ti.end
        assert loaded.all_starts() == ti.all_starts()
        assert loaded.all_targets() == ti.all_targets()
        np.testing.assert_allclose(
            [h for hs in loaded.all_tmrcas() for h in hs],
            [h for hs in ti.all_tmrcas() for h in hs],
        )
        assert loaded.all_mismatches() == ti.all_mismatches()

    def test_serialize_with_metadata(self):
        import pandas as pd
        from threads_arg import ThreadingInstructions
        from threads_arg.serialization import serialize_instructions, load_instructions, load_metadata, load_sample_names

        positions = [100, 200, 300]
        # 2 haploid samples (1 diploid) so sample_names has 1 entry
        ti = ThreadingInstructions(
            [[0], [0]], [[50.0], [60.0]], [[1], [0]], [[], []], positions, 0, 400
        )

        metadata = pd.DataFrame({
            "CHROM": ["1", "1", "1"],
            "POS": [100, 200, 300],
            "ID": ["rs1", "rs2", "rs3"],
            "REF": ["A", "C", "G"],
            "ALT": ["T", "G", "A"],
            "QUAL": [".", ".", "."],
            "FILTER": ["PASS", "PASS", "PASS"],
        })
        sample_names = ["sample_0"]

        with tempfile.NamedTemporaryFile(suffix=".threads", delete=False) as f:
            out_path = f.name

        serialize_instructions(ti, out_path, variant_metadata=metadata, sample_names=sample_names)

        loaded_meta = load_metadata(out_path)
        assert list(loaded_meta.columns) == ["CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER"]
        assert len(loaded_meta) == 3

        loaded_names = load_sample_names(out_path)
        # load_sample_names returns bytes from HDF5
        assert [n.decode() if isinstance(n, bytes) else n for n in loaded_names] == ["sample_0"]
