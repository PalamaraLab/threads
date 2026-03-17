"""
Tests for the `threads vcf` command: GenotypeIterator, VCFWriter, and threads_to_vcf.
"""
import ctypes
import os
import re
import sys
import tempfile

import numpy as np
import pytest

# Flush all C stdio streams (needed to capture C++ stdout from VCFWriter)
_libc = ctypes.CDLL(None)
def _flush_c_stdout():
    _libc.fflush(None)

from pathlib import Path

# ---------------------------------------------------------------------------
# Paths
# ---------------------------------------------------------------------------
TEST_DATA = Path(__file__).parent / "data"
THREADS_FIT = TEST_DATA / "expected_infer_fit_to_data_snapshot.threads"
THREADS_NO_FIT = TEST_DATA / "expected_infer_snapshot.threads"
PANEL_PGEN = TEST_DATA / "panel.pgen"
PANEL_PVAR = TEST_DATA / "panel.pvar"
PANEL_PSAM = TEST_DATA / "panel.psam"


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------
def _load_instructions(threads_path):
    from threads_arg.serialization import load_instructions
    return load_instructions(str(threads_path))


def _load_sample_names():
    names = []
    with open(PANEL_PSAM) as f:
        for line in f:
            if line.startswith("#"):
                continue
            names.append(line.strip().split()[0])
    return names


def _capture_vcf_output(instructions, variant_metadata, sample_names):
    """Run VCFWriter.write_vcf() and capture C++ stdout via fd redirect."""
    from threads_arg import VCFWriter

    with tempfile.NamedTemporaryFile(mode="w+", suffix=".vcf", delete=False) as tmp:
        tmpname = tmp.name

    try:
        old_fd = os.dup(1)
        new_fd = os.open(tmpname, os.O_WRONLY | os.O_TRUNC | os.O_CREAT)
        os.dup2(new_fd, 1)
        os.close(new_fd)

        writer = VCFWriter(instructions)
        writer.set_chrom(variant_metadata["CHROM"].astype(str))
        writer.set_pos(variant_metadata["POS"].astype(str))
        writer.set_id(variant_metadata["ID"].astype(str))
        writer.set_ref(variant_metadata["REF"].astype(str))
        writer.set_alt(variant_metadata["ALT"].astype(str))
        writer.set_qual(variant_metadata["QUAL"].astype(str))
        writer.set_filter(variant_metadata["FILTER"].astype(str))
        writer.set_sample_names(sample_names)
        writer.write_vcf()

        _flush_c_stdout()
        sys.stdout.flush()
        os.dup2(old_fd, 1)
        os.close(old_fd)

        with open(tmpname) as f:
            return f.readlines()
    finally:
        os.unlink(tmpname)


@pytest.fixture(scope="module")
def fit_instructions():
    return _load_instructions(THREADS_FIT)


@pytest.fixture(scope="module")
def nofit_instructions():
    return _load_instructions(THREADS_NO_FIT)


@pytest.fixture(scope="module")
def variant_metadata():
    from threads_arg.utils import read_variant_metadata
    return read_variant_metadata(str(PANEL_PGEN))


@pytest.fixture(scope="module")
def sample_names():
    return _load_sample_names()


@pytest.fixture(scope="module")
def panel_genotypes():
    """All genotypes from panel.pgen as (n_sites, n_haps) int32 array."""
    from threads_arg.utils import read_all_genotypes
    return read_all_genotypes(str(PANEL_PGEN))


@pytest.fixture(scope="module")
def vcf_lines(fit_instructions, variant_metadata, sample_names):
    return _capture_vcf_output(fit_instructions, variant_metadata, sample_names)


# ===================================================================
# GenotypeIterator unit tests
# ===================================================================
class TestGenotypeIterator:

    def test_genotype_values_biallelic(self, fit_instructions):
        """All genotype values should be 0 or 1."""
        from threads_arg import GenotypeIterator
        gi = GenotypeIterator(fit_instructions)
        for _ in range(50):
            g = gi.next_genotype()
            assert set(g).issubset({0, 1}), f"non-biallelic values: {set(g) - {0, 1}}"

    def test_genotype_length_matches_num_samples(self, fit_instructions):
        """Each genotype vector should have length == num_samples (haploid)."""
        from threads_arg import GenotypeIterator
        gi = GenotypeIterator(fit_instructions)
        g = gi.next_genotype()
        assert len(g) == fit_instructions.num_samples

    def test_total_sites(self, fit_instructions):
        """Iterator should yield exactly num_sites genotypes."""
        from threads_arg import GenotypeIterator
        gi = GenotypeIterator(fit_instructions)
        count = 0
        while gi.has_next_genotype():
            gi.next_genotype()
            count += 1
        assert count == fit_instructions.num_sites

    def test_round_trip_fit_to_data(self, fit_instructions, panel_genotypes):
        """fit_to_data .threads should exactly reconstruct panel.pgen genotypes."""
        from threads_arg import GenotypeIterator
        gi = GenotypeIterator(fit_instructions)
        for site_idx in range(fit_instructions.num_sites):
            g = np.array(gi.next_genotype())
            np.testing.assert_array_equal(
                g, panel_genotypes[site_idx],
                err_msg=f"mismatch at site {site_idx}"
            )

    def test_round_trip_no_fit(self, nofit_instructions, panel_genotypes):
        """non-fit_to_data .threads: check genotype reconstruction."""
        from threads_arg import GenotypeIterator
        gi = GenotypeIterator(nofit_instructions)
        mismatches = 0
        for site_idx in range(nofit_instructions.num_sites):
            g = np.array(gi.next_genotype())
            if not np.array_equal(g, panel_genotypes[site_idx]):
                mismatches += 1
        # This snapshot also has zero mismatches empirically
        assert mismatches == 0, f"{mismatches}/{nofit_instructions.num_sites} sites differ"

    def test_allele_frequency_reasonable(self, fit_instructions, panel_genotypes):
        """Allele frequencies from iterator should match panel.pgen."""
        from threads_arg import GenotypeIterator
        gi = GenotypeIterator(fit_instructions)
        for site_idx in range(min(100, fit_instructions.num_sites)):
            g = np.array(gi.next_genotype())
            assert np.mean(g) == pytest.approx(
                np.mean(panel_genotypes[site_idx]), abs=1e-10
            )


# ===================================================================
# VCFWriter format tests
# ===================================================================
class TestVCFWriterFormat:

    def test_header_starts_vcf42(self, vcf_lines):
        assert vcf_lines[0].strip() == "##fileformat=VCFv4.2"

    def test_has_source_line(self, vcf_lines):
        assert any(l.startswith("##source=") for l in vcf_lines)

    def test_has_contig_line(self, vcf_lines):
        assert any(l.startswith("##contig=") for l in vcf_lines)

    def test_chrom_header_present(self, vcf_lines):
        header_lines = [l for l in vcf_lines if l.startswith("#CHROM")]
        assert len(header_lines) == 1

    def test_chrom_header_columns(self, vcf_lines, sample_names):
        for l in vcf_lines:
            if l.startswith("#CHROM"):
                fields = l.strip().split("\t")
                fixed = fields[:9]
                assert fixed == [
                    "#CHROM", "POS", "ID", "REF", "ALT",
                    "QUAL", "FILTER", "INFO", "FORMAT"
                ]
                # Sample columns
                sample_cols = fields[9:]
                assert len(sample_cols) == len(sample_names)
                assert sample_cols[0] == sample_names[0]
                assert sample_cols[-1] == sample_names[-1]
                break

    def test_data_line_count(self, vcf_lines, fit_instructions):
        """Number of data lines should equal number of sites."""
        data_lines = [l for l in vcf_lines if not l.startswith("#")]
        assert len(data_lines) == fit_instructions.num_sites

    def test_format_field_is_gt(self, vcf_lines):
        """FORMAT column should be 'GT'."""
        for l in vcf_lines:
            if not l.startswith("#"):
                fields = l.split("\t")
                assert fields[8] == "GT"
                break

    def test_genotype_format_phased(self, vcf_lines):
        """Genotype fields should match pattern [01]|[01]."""
        gt_pattern = re.compile(r"^[01]\|[01]$")
        for l in vcf_lines:
            if not l.startswith("#"):
                fields = l.rstrip("\n").split("\t")
                # Skip empty trailing field from known trailing-tab bug
                gt_fields = [f for f in fields[9:] if f]
                for gt in gt_fields[:20]:  # spot-check first 20
                    assert gt_pattern.match(gt), f"bad GT format: {gt!r}"
                break

    def test_info_field_has_ns(self, vcf_lines, sample_names):
        """INFO field should contain NS=<num_diploid_samples>."""
        expected_ns = f"NS={len(sample_names)}"
        for l in vcf_lines:
            if not l.startswith("#"):
                fields = l.split("\t")
                assert fields[7] == expected_ns
                break

    def test_trailing_tab_known_bug(self, vcf_lines):
        """Document: VCFWriter emits a trailing tab on each data line."""
        for l in vcf_lines:
            if not l.startswith("#"):
                # The raw line (before rstrip) should end with \t\n
                assert l.endswith("\t\n") or l.endswith("\t"), \
                    "trailing tab bug may have been fixed — update this test"
                break


# ===================================================================
# threads_to_vcf integration tests
# ===================================================================
class TestThreadsToVcf:

    def test_missing_metadata_raises(self):
        """Calling threads_to_vcf without external files on a .threads
        that lacks embedded metadata should raise RuntimeError."""
        from threads_arg.threads_to_vcf import threads_to_vcf
        with pytest.raises(RuntimeError, match="Unable to load sample information"):
            threads_to_vcf(str(THREADS_FIT))

    def test_with_external_pvar(self, fit_instructions, sample_names):
        """threads_to_vcf with external .pvar and samples file should produce output."""
        from threads_arg.threads_to_vcf import threads_to_vcf

        with tempfile.NamedTemporaryFile(mode="w", suffix=".txt", delete=False) as sf:
            for name in sample_names:
                sf.write(name + "\n")
            samples_path = sf.name

        with tempfile.NamedTemporaryFile(mode="w+", suffix=".vcf", delete=False) as tmp:
            tmpname = tmp.name

        try:
            old_fd = os.dup(1)
            new_fd = os.open(tmpname, os.O_WRONLY | os.O_TRUNC | os.O_CREAT)
            os.dup2(new_fd, 1)
            os.close(new_fd)

            threads_to_vcf(
                str(THREADS_FIT),
                samples=samples_path,
                variants=str(PANEL_PVAR),
            )

            _flush_c_stdout()
            sys.stdout.flush()
            os.dup2(old_fd, 1)
            os.close(old_fd)

            with open(tmpname) as f:
                lines = f.readlines()

            assert len(lines) > 0, "no output produced"
            data_lines = [l for l in lines if not l.startswith("#")]
            assert len(data_lines) == fit_instructions.num_sites
        finally:
            os.unlink(tmpname)
            os.unlink(samples_path)

    def test_genotype_consistency(self, variant_metadata, sample_names, panel_genotypes):
        """Genotypes in VCF output should match panel.pgen."""
        instructions = _load_instructions(THREADS_FIT)
        lines = _capture_vcf_output(instructions, variant_metadata, sample_names)

        n_diploid = len(sample_names)
        # Check first 100 data lines
        site_idx = 0
        for l in lines:
            if l.startswith("#"):
                continue
            if site_idx >= 100:
                break
            fields = l.rstrip("\n").split("\t")
            gt_fields = [f for f in fields[9:] if f]  # skip trailing empty
            assert len(gt_fields) == n_diploid

            # Reconstruct haploid genotypes from phased GTs
            haps = np.zeros(2 * n_diploid, dtype=int)
            for i, gt_str in enumerate(gt_fields):
                a, b = gt_str.split("|")
                haps[2 * i] = int(a)
                haps[2 * i + 1] = int(b)

            np.testing.assert_array_equal(
                haps, panel_genotypes[site_idx],
                err_msg=f"VCF genotype mismatch at site {site_idx}"
            )
            site_idx += 1
