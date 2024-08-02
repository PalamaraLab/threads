# Threads

Highly scalable inference of ancestral recombination graphs (ARGs)

## Installation

```sh
pip install threads_arg
```

## Usage

### ARG inference

You will need

- genotypes in pgen format
- list of variants in bim or pvar format (with the same prefix as the pgen)
- genetic map with 4 columns: Chromosome, SNP, cM, bp
- demography file with two columns: generations in the past, effective population size in haploids

Minimal usage using the provided example data:

```sh
threads infer \
    --pgen example/example_data.pgen \
    --map_gz example/example_data.map \
    --demography example/Ne10000.demo \
    --out example/example_data.threads

threads convert \
    --threads example/example_data.threads \
    --argn example/example_data.argn
```

This will write a `.threads` file to `path/to/output.threads`.

`threads infer` accepts more options:

```sh
threads infer \
    --pgen path/to/input.pgen \
    --map_gz path/to/genetic_map.gz \
    --demography path/to/demography \
    --out path/to/output.threads \
    --modality [wgs|array] (default: wgs) \
    --query_interval (default: 0.01) \
    --match_group_interval (default: 0.5) \
    --max_sample_batch_size (default: None) \
    --mutation_rate (default: 1.4e-8) \
    --region 1234-56789 (default: whole region, end-inclusive) \
    --num_threads 8 (default: 1)
```

`--modality array` can be set for inference from arrays. 

`--query_interval` and `--match_group_interval` can be raised to save memory for inference over long genomic regions, this will have little impact on accuracy, especially for sparse variants. 

The HMM mutation rate can be set with `--mutation_rate`. This defaults to a realistic human rate of `1.4e-8` per site per generation.

Specifying a `--region start-end` means the output ARG is truncated to those base-pair coordinates (end-inclusive). The whole input set will still be used for inference.

Parallelism can be enabled by specifying `--num_threads`

### ARG conversion

`.threads` files can be converted to `.argn` and `.tsz` using

```sh
threads convert \
    --threads arg.threads \
    --argn arg.argn
```

and

```sh
threads convert \
    --threads arg.threads \
    --tsz arg.tsz
```

### Phasing/imputation/variant mapping

These functions are in an experimental stage and will be released later.
