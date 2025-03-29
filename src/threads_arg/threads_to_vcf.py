from .serialization import load_instructions, load_metadata, load_sample_names
from .utils import read_variant_metadata
from threads_arg import VCFWriter

def threads_to_vcf(threads, samples=None, variants=None):
    if samples is None:
        try:
            sample_names = load_sample_names(threads)
        except KeyError:
            raise RuntimeError("Unable to load sample information from threading instructions. This may because files were inferred using an older version of threads. Please provide files with variant information (in .bim/.pvar format) and sample IDs in a file with one sample per line.")
    else:
        with open(samples, "r") as samplefile:
            sample_names = [l.strip() for l in samplefile]

    if variants is None:
        try:
            variant_metadata = load_metadata(threads)
        except KeyError:
            raise RuntimeError("Unable to load sample information from threading instructions. This may because files were inferred using an older version of threads. Please provide files with variant information (in .bim/.pvar format) and sample IDs in a file with one sample per line.")
    else:
        assert variants.endswith(".bim") or variants.endswith(".var")
        if variants.endswith(".bim"):
            variant_metadata = read_variant_metadata(variants[:-4] + ".pgen")
        if variants.endswith(".pvar"):
            variant_metadata = read_variant_metadata(variants[:-5] + ".pgen")

    instructions = load_instructions(threads)
    assert len(sample_names) == len(instructions.all_starts()) // 2
    assert variant_metadata.shape[0] == len(instructions.positions)

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
