"""
Tests for the `threads map` command: _mapping_string formatting, get_leaf_ids_at,
MAF filtering logic, and full pipeline validation.
"""
import tempfile

import numpy as np
import pytest
import arg_needle_lib

from pathlib import Path

TEST_DATA = Path(__file__).parent / "data"
ARGN_SNAPSHOT = str(TEST_DATA / "expected_convert_snapshot.argn")
MUT_SNAPSHOT = str(TEST_DATA / "expected_mapping_snapshot.mut")
PANEL_VCF = str(TEST_DATA / "panel.vcf.gz")
REGION = "1:400000-600000"


# ===================================================================
# _mapping_string: format mutation-to-edge mappings
# ===================================================================
class TestMappingString:
    @pytest.fixture()
    def mock_edge(self):
        """Create a minimal mock edge with child/parent height attributes."""
        class MockNode:
            def __init__(self, height):
                self.height = height
                self.ID = 0
        class MockEdge:
            def __init__(self, child_h, parent_h):
                self.child = MockNode(child_h)
                self.parent = MockNode(parent_h)
        return MockEdge

    def test_empty_edges_returns_nan(self):
        from threads_arg.map_mutations_to_arg import _mapping_string
        assert _mapping_string([], []) == "NaN"

    def test_single_edge_format(self, mock_edge):
        from threads_arg.map_mutations_to_arg import _mapping_string
        edge = mock_edge(0.0, 1138.7262)
        result = _mapping_string([[-1]], [edge])
        assert result == "-1,0.0000,1138.7262"

    def test_single_edge_precision(self, mock_edge):
        from threads_arg.map_mutations_to_arg import _mapping_string
        edge = mock_edge(378.9087, 1303.8883)
        result = _mapping_string([[-1]], [edge])
        assert result == "-1,378.9087,1303.8883"

    def test_multiple_edges_format(self, mock_edge):
        from threads_arg.map_mutations_to_arg import _mapping_string
        edges = [mock_edge(0.0, 599.4252), mock_edge(0.0, 659.2068), mock_edge(0.0, 704.7969)]
        carrier_sets = [[124], [740], [816]]
        result = _mapping_string(carrier_sets, edges)
        parts = result.split(";")
        assert len(parts) == 3
        assert parts[0] == "124,0.0000,599.4252"
        assert parts[1] == "740,0.0000,659.2068"
        assert parts[2] == "816,0.0000,704.7969"

    def test_single_edge_always_uses_sentinel(self, mock_edge):
        """A single edge always uses -1 sentinel regardless of carrier set."""
        from threads_arg.map_mutations_to_arg import _mapping_string
        edges = [mock_edge(100.0, 500.0)]
        carrier_sets = [[5, 10, 15]]
        result = _mapping_string(carrier_sets, edges)
        assert result == "-1,100.0000,500.0000"

    def test_two_edges_with_mixed_carriers(self, mock_edge):
        from threads_arg.map_mutations_to_arg import _mapping_string
        edges = [mock_edge(693.77, 918.05), mock_edge(0.0, 1618.53)]
        carrier_sets = [[569, 147, 72], [156]]
        result = _mapping_string(carrier_sets, edges)
        parts = result.split(";")
        assert len(parts) == 2
        assert parts[0].startswith("569.147.72,")
        assert parts[1].startswith("156,")


# ===================================================================
# get_leaf_ids_at: recursive leaf extraction from ARG edges
# ===================================================================
class TestGetLeafIds:
    @pytest.fixture(scope="class")
    def arg_with_roots(self):
        arg = arg_needle_lib.deserialize_arg(ARGN_SNAPSHOT)
        arg.populate_children_and_roots()
        return arg

    def test_leaf_ids_non_empty(self, arg_with_roots):
        """Mapping a real variant should produce non-empty leaf IDs."""
        from threads_arg.map_mutations_to_arg import get_leaf_ids_at
        from cyvcf2 import VCF

        vcf = VCF(PANEL_VCF)
        for record in vcf(REGION):
            ac = int(record.INFO.get("AC"))
            an = int(record.INFO.get("AN"))
            maf = min(ac / an, 1 - ac / an)
            if maf == 0 or maf > 0.01:
                continue

            hap = np.array(record.genotypes)[:, :2].flatten()
            mapping, _ = arg_needle_lib.map_genotype_to_ARG_approximate(
                arg_with_roots, hap, float(record.POS - arg_with_roots.offset)
            )
            if len(mapping) > 1:
                for edge in mapping:
                    leaves = get_leaf_ids_at(arg_with_roots, edge, record.POS)
                    assert len(leaves) > 0, "edge should subtend at least one leaf"
                return  # tested one multi-edge mapping, done
        pytest.skip("No multi-edge mapping found in region")

    def test_leaf_ids_are_leaves(self, arg_with_roots):
        """All returned IDs should be actual leaf nodes."""
        from threads_arg.map_mutations_to_arg import get_leaf_ids_at
        from cyvcf2 import VCF

        vcf = VCF(PANEL_VCF)
        for record in vcf(REGION):
            ac = int(record.INFO.get("AC"))
            an = int(record.INFO.get("AN"))
            maf = min(ac / an, 1 - ac / an)
            if maf == 0 or maf > 0.01:
                continue

            hap = np.array(record.genotypes)[:, :2].flatten()
            mapping, _ = arg_needle_lib.map_genotype_to_ARG_approximate(
                arg_with_roots, hap, float(record.POS - arg_with_roots.offset)
            )
            if len(mapping) > 1:
                for edge in mapping:
                    leaves = get_leaf_ids_at(arg_with_roots, edge, record.POS)
                    for lid in leaves:
                        assert arg_with_roots.is_leaf(lid)
                return
        pytest.skip("No multi-edge mapping found in region")


# ===================================================================
# .mut file format validation
# ===================================================================
class TestMutFileFormat:
    @pytest.fixture(scope="class")
    def mut_lines(self):
        with open(MUT_SNAPSHOT) as f:
            return f.readlines()

    def test_non_empty(self, mut_lines):
        assert len(mut_lines) > 0

    def test_tab_separated_four_columns(self, mut_lines):
        for i, line in enumerate(mut_lines[:50]):
            fields = line.rstrip("\n").split("\t")
            assert len(fields) == 4, f"line {i+1}: expected 4 tab-separated fields, got {len(fields)}"

    def test_position_column_is_int(self, mut_lines):
        for line in mut_lines:
            fields = line.rstrip("\n").split("\t")
            int(fields[1])  # should not raise

    def test_flipped_column_is_binary(self, mut_lines):
        for line in mut_lines:
            fields = line.rstrip("\n").split("\t")
            assert fields[2] in ("0", "1")

    def test_mapping_string_format(self, mut_lines):
        """Mapping string should be NaN, or -1,h1,h2, or semicolon-separated groups."""
        for line in mut_lines:
            fields = line.rstrip("\n").split("\t")
            ms = fields[3]
            if ms == "NaN":
                continue
            groups = ms.split(";")
            for group in groups:
                parts = group.split(",")
                assert len(parts) == 3, f"bad mapping group: {group!r}"
                # Last two should be floats
                float(parts[1])
                float(parts[2])

    def test_positions_in_region(self, mut_lines):
        """All positions should be within the expected region."""
        for line in mut_lines:
            fields = line.rstrip("\n").split("\t")
            pos = int(fields[1])
            assert 400000 <= pos <= 600000

    def test_uniquely_mapped_uses_sentinel(self, mut_lines):
        """Single-edge mappings should use -1 sentinel."""
        for line in mut_lines:
            fields = line.rstrip("\n").split("\t")
            ms = fields[3]
            if ms == "NaN":
                continue
            if ";" not in ms:
                # Single edge
                assert ms.startswith("-1,"), f"single-edge mapping should start with -1: {ms!r}"


# ===================================================================
# Full pipeline: threads_map_mutations_to_arg
# ===================================================================
class TestMapPipeline:
    def test_produces_mut_file(self):
        from threads_arg.map_mutations_to_arg import threads_map_mutations_to_arg
        with tempfile.TemporaryDirectory() as tmpdir:
            out_mut = str(Path(tmpdir) / "test.mut")
            threads_map_mutations_to_arg(
                argn=ARGN_SNAPSHOT, out=out_mut, maf=0.01,
                input=PANEL_VCF, region=REGION, num_threads=1
            )
            assert Path(out_mut).exists()
            with open(out_mut) as f:
                lines = f.readlines()
            assert len(lines) > 0

    def test_same_variants_mapped(self):
        """Same variant IDs and positions should be mapped as the snapshot."""
        from threads_arg.map_mutations_to_arg import threads_map_mutations_to_arg
        with tempfile.TemporaryDirectory() as tmpdir:
            out_mut = str(Path(tmpdir) / "test.mut")
            threads_map_mutations_to_arg(
                argn=ARGN_SNAPSHOT, out=out_mut, maf=0.01,
                input=PANEL_VCF, region=REGION, num_threads=1
            )
            with open(out_mut) as gen, open(MUT_SNAPSHOT) as exp:
                gen_lines = gen.readlines()
                exp_lines = exp.readlines()

            assert len(gen_lines) == len(exp_lines), \
                f"line count: {len(gen_lines)} vs {len(exp_lines)}"
            # Variant IDs and positions should match exactly
            for i, (gl, el) in enumerate(zip(gen_lines, exp_lines)):
                g_fields = gl.split("\t")
                e_fields = el.split("\t")
                assert g_fields[0] == e_fields[0], f"line {i+1}: variant ID differs"
                assert g_fields[1] == e_fields[1], f"line {i+1}: position differs"
                assert g_fields[2] == e_fields[2], f"line {i+1}: flipped flag differs"

    def test_stricter_maf_fewer_variants(self):
        """Tighter MAF filter should produce fewer or equal mappings."""
        from threads_arg.map_mutations_to_arg import threads_map_mutations_to_arg
        with tempfile.TemporaryDirectory() as tmpdir:
            out_wide = str(Path(tmpdir) / "wide.mut")
            out_strict = str(Path(tmpdir) / "strict.mut")

            threads_map_mutations_to_arg(
                argn=ARGN_SNAPSHOT, out=out_wide, maf=0.05,
                input=PANEL_VCF, region=REGION, num_threads=1
            )
            threads_map_mutations_to_arg(
                argn=ARGN_SNAPSHOT, out=out_strict, maf=0.005,
                input=PANEL_VCF, region=REGION, num_threads=1
            )
            with open(out_wide) as f:
                n_wide = len(f.readlines())
            with open(out_strict) as f:
                n_strict = len(f.readlines())
            assert n_strict <= n_wide
