# Installation (rescomp only):

Load everything we'll need
```
module load GCC/11.3.0
module load CMake/3.23.1-GCCcore-11.3.0
module load git/2.36.0-GCCcore-11.3.0-nodocs
module load Python/3.10.4-GCCcore-11.3.0
module load GSL/2.7-GCC-11.3.0
module load Boost/1.79.0-GCC-11.3.0
module load Eigen/3.4.0-GCCcore-11.3.0
module load pybind11/2.9.2-GCCcore-11.3.0
```

NB: Eigen is only used for functionality in TGEN.cpp. We can consider releasing without TGEN to have fewer dependencies.

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
```
# Run the inference
threads infer \
    --pgen path/to/input.pgen \
    --map_gz path/to/genetic_map.gz \
    --demography path/to/demography \
    --out path/to/output.threads \
    --modality [wgs|array] (default: wgs) \
    --mutation_rate (default: 1.4e-8) \
    --region 1234-56789 (default: whole region, end-inclusive) \
    --num_threads 8 (default: 1)
```
This will write a `.threads` file to `path/to/output.threads`.

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
