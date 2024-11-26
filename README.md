# Installation (rescomp only):

Load everything we'll need
```
module load GCC/11.3.0
module load CMake/3.23.1-GCCcore-11.3.0
module load git/2.36.0-GCCcore-11.3.0-nodocs
module load Python/3.10.4-GCCcore-11.3.0
module load GSL/2.7-GCC-11.3.0
module load Boost/1.79.0-GCC-11.3.0
module load pybind11/2.9.2-GCCcore-11.3.0
```

Fire up and activate a new venv:
```
python -m venv venv
. venv/bin/activate

pip install --upgrade pip setuptools wheel
pip install cmake ninja
```

Then install Threads:
```
git clone https://github.com/PalamaraLab/TDPBWT.git
cd TDPBWT
pip install .
```

# Usage
## ARG inference

You will need
- genotypes in pgen format
- list of variants in bim or pvar format (with the same prefix as the pgen)
- genetic map with 4 columns: Chromosome, SNP, cM, bp
- demography file with two columns: generations in the past, effective population size in haploids

Minimal usage using the provided example data:
```
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
```
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
    --num_threads 8 (default: 1) \
    --fit_to_data (default: False)
    --allele_ages (default: None)
```

`--modality array` can be set for inference from arrays. 

`--query_interval` and `--match_group_interval` can be raised to save memory for inference over long genomic regions, this will have little impact on accuracy, especially for sparse variants. 

The HMM mutation rate can be set with `--mutation_rate`. This defaults to a realistic human rate of `1.4e-8` per site per generation.

Specifying a `--region start-end` means the output ARG is truncated to those base-pair coordinates (end-inclusive). The whole input set will still be used for inference.

Parallelism can be enabled by specifying `--num_threads`.

By specifying `--fit_to_data`, the inferred threading instructions are manipulated post-hoc to ensure that each input mutation is uniquely represented by an edge in the ARG. This may lead to more fragmented ARGs.

If `--fit_to_data` is specified and a file with two columns, a SNP ID and an allele age, is provided with `--allele_ages`, the ARG is altered in a way that makes it consistent with those mutations. If `--allele_ages` is not specified, allele ages are automatically inferred. If `--fit_to_data` is not specified, this option is ignored.

## ARG conversion
`.threads` files can be converted to `.argn` and `.tsz` using
```
threads convert \
    --threads arg.threads \
    --argn arg.argn
```
and 
```
threads convert \
    --threads arg.threads \
    --tsz arg.tsz
```

## Phasing/imputation/variant mapping
These functions are in an experimental phase and will be released later.

## Allele age estimation
Threads can be used to infer the age of any set of genotypes overlapping the genomic region covered by the ARG using the command `threads allele-ages`:
```
threads allele-ages --threads input.threads \
    --pgen input.pgen \
    --region 1:234-567 (Optional) \
    --out output_4/allele.ages
```

## Data consistency
The manipulation of threading instructions performed using the `--fit_to_data` flag in `threads infer` can be performed separately using the command `threads fit-to-data`:
```
threads fit-to-data --threads input.threads \
    --pgen input.pgen \
    --allele_ages input.ages (Optional) \
    --region 1:234-567 (Optional) \
    --out output_4/allele.ages
```
The threading instructions can be made to fit to any genotypes provided with the `input.pgen`.