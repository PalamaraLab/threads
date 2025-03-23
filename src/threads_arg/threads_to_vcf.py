from .serialization import load_instructions, load_metadata, load_sample_names
from threads_arg import VCFWriter

def threads_to_vcf(threads):
    instructions = load_instructions(threads)
    sample_names = load_sample_names(threads)
    variant_metadata = load_metadata(threads)
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
