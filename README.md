Steps towards set-up:

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
module load htslib
# for boost_iostreams
LD_LIBRARY_PATH=/users/palamara/awo066/lib:$LD_LIBRARY_PATH
```

Fire up and activate a new venv:
```
python -m venv venv
. venv/bin/activate

pip install --upgrade pip setuptools wheel
pip install cmake ninja
```

`arg_needle_lib` is not yet on pypi so we pull and install it:
```
git clone https://github.com/PalamaraLab/arg_needle_lib.git
cd arg_needle_lib
pip install .
cd ..
```

Then do the same for threads:
```
git clone https://github.com/PalamaraLab/TDPBWT.git
cd TDPBWT
pip install .
```

Make sure it works:
```
cd examples/example_data
# Convert genotypes into zarr archives
threads files --vcf data.[vcf/vcf.gz/bcf] --zarr data.zarr
# Run the inference
threads infer \
    --genotypes snp_data \
    --map_gz snp_data.map.gz \
    --mutation_rate 1.65e-8 \
    --demography CEU.demo \
    --modality [array/seq] \
    --threads snp_data.threads
# Convert to .argn and/or .tsz
threads convert \
    --threads snp_data.threads \
    --argn snp_data.argn \
    --tsz snp_data.tsz

# And for sequencing data:
threads files \
    --hap_gz seq_data.hap.gz \
    --sample seq_data.sample \
    --bfile seq_data
threads infer \
    --bfile seq_data \
    --map_gz seq_data.map.gz \
    --mutation_rate 1.65e-8 \
    --demography CEU.demo \
    --modality seq \
    --threads seq_data.threads
threads convert \
    --threads seq_data.threads \
    --argn seq_data.argn \
    --tsz seq_data.tsz
