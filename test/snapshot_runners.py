from pathlib import Path


BASE_DIR = Path(__file__).parent.parent
TEST_DATA_DIR = BASE_DIR / "test" / "data"


def run_infer_snapshot(out_threads_path: Path):
    from threads_arg.infer import threads_infer
    threads_infer(
        pgen=str(TEST_DATA_DIR / "panel.pgen"),
        map=str(TEST_DATA_DIR / "gmap_02.map"),
        recombination_rate=1.3e-8,
        demography=str(TEST_DATA_DIR / "CEU_unscaled.demo"),
        mutation_rate=1.4e-8,
        query_interval=0.01,
        match_group_interval=0.5,
        mode="wgs",
        num_threads=1,
        region=None,
        fit_to_data=False,
        allele_ages=None,
        max_sample_batch_size=None,
        out=str(out_threads_path)
    )


def run_convert_snapshot(in_threads_path: Path, out_argn_path: Path):
    from threads_arg.convert import threads_convert
    threads_convert(
        threads=str(in_threads_path),
        argn=str(out_argn_path),
        tsz=None
    )


def run_map_snapshot(in_argn_path: Path, out_mut_path: Path):
    from threads_arg.map_mutations_to_arg import threads_map_mutations_to_arg
    threads_map_mutations_to_arg(
        argn=str(in_argn_path),
        out=str(out_mut_path),
        maf=0.01,
        input=str(TEST_DATA_DIR / "panel.vcf.gz"),
        region="1:400000-600000",
        num_threads=1
    )


def run_impute_snapshot(in_mut_path: Path, out_vcf_path: Path):
    from threads_arg.impute import Impute
    Impute(
        TEST_DATA_DIR / "panel.vcf.gz",
        TEST_DATA_DIR / "target.vcf.gz",
        TEST_DATA_DIR / "gmap_04.map",
        in_mut_path,
        TEST_DATA_DIR / "CEU_unscaled.demo",
        out_vcf_path,
        "1:400000-600000"
    )
