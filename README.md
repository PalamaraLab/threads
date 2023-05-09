Installation (rescomp only):

Load everything we'll need (some of these won't be necessary but I haven't checked which ones) 
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

Fire up and activate a new venv:
```
python -m venv venv
. venv/bin/activate

pip install --upgrade pip setuptools wheel
pip install cmake ninja
```

Then do the same for threads:
```
git clone https://github.com/PalamaraLab/TDPBWT.git
cd TDPBWT
pip install .
```

To make sure it works you need
- genotypes in pgen format
- genetic map with 4 columns: Chromosome, SNP_id, cM, bp and one row for each marker in the pgen
- demography file, interpreted as haploid, with two columns: generation, effective population size
```
# Run the inference
threads infer \
    --pgen path/to/input.pgen \
    --map_gz path/to/genetic_map.gz \
    --mutation_rate 1.65e-8 \
    --demography path/to/demography \
    --modality [array|wgs] \
    --out path/to/output.threads
    --region 123-456789 (optional)
    --num_threads 8 (optional)

# Convert to .argn and/or .tsz
threads convert \
    --threads snp_data.threads \
    --argn snp_data.argn \
    --tsz snp_data.tsz
```