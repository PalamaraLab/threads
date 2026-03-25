from .serialization import load_instructions, load_metadata, load_sample_names
from .utils import read_variant_metadata
import numpy as np
import pgenlib


def threads_to_pgen(threads, out, samples=None, variants=None):
    """Write threading instructions to pgen/pvar/psam files.

    Args:
        threads: Path to threading instructions file.
        out: Output prefix (will create out.pgen, out.pvar, out.psam).
        samples: Optional path to sample names file (one per line).
        variants: Optional path to .bim or .pvar with variant metadata.
    """
    if samples is None:
        try:
            sample_names = load_sample_names(threads)
        except (KeyError, RuntimeError):
            raise RuntimeError(
                "Unable to load sample information. Provide a samples file.")
    else:
        with open(samples) as f:
            sample_names = [line.strip() for line in f]

    if variants is None:
        try:
            variant_metadata = load_metadata(threads)
        except (KeyError, RuntimeError):
            raise RuntimeError(
                "Unable to load variant metadata. Provide a .bim or .pvar file.")
    else:
        assert variants.endswith(".bim") or variants.endswith(".pvar")
        stem = variants.rsplit(".", 1)[0]
        variant_metadata = read_variant_metadata(stem + ".pgen")

    instructions = load_instructions(threads)
    n_haps = instructions.num_samples
    n_dip = n_haps // 2
    n_sites = instructions.num_sites

    if len(sample_names) != n_dip:
        raise RuntimeError(
            f"Sample count mismatch: {len(sample_names)} names vs "
            f"{n_dip} diploid samples in instructions.")
    if variant_metadata.shape[0] != n_sites:
        raise RuntimeError(
            f"Variant count mismatch: {variant_metadata.shape[0]} metadata rows vs "
            f"{n_sites} sites in instructions.")

    # Write pgen
    geno = instructions.genotype_matrix_numpy()  # (n_sites, n_haps) int8
    with pgenlib.PgenWriter(
            filename=(out + ".pgen").encode(),
            sample_ct=n_dip,
            variant_ct=n_sites,
            hardcall_phase_present=True) as writer:
        alleles = np.empty(n_haps, dtype=np.int32)
        for site in range(n_sites):
            alleles[:] = geno[site]
            writer.append_alleles(alleles, all_phased=True)

    # Write pvar
    with open(out + ".pvar", "w") as f:
        f.write("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\n")
        chroms = variant_metadata["CHROM"]
        positions = variant_metadata["POS"]
        ids = variant_metadata["ID"]
        refs = variant_metadata["REF"]
        alts = variant_metadata["ALT"]
        quals = variant_metadata["QUAL"]
        filters = variant_metadata["FILTER"]
        for i in range(n_sites):
            f.write(f"{chroms[i]}\t{positions[i]}\t{ids[i]}\t"
                    f"{refs[i]}\t{alts[i]}\t{quals[i]}\t{filters[i]}\n")

    # Write psam
    with open(out + ".psam", "w") as f:
        f.write("#IID\n")
        for name in sample_names:
            f.write(f"{name}\n")
