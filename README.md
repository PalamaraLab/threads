Steps towards set-up:

Load everything we'll need (some of these won't be necessary but I haven't checked which ones) 
```
module load GCC/10.2.0
module load CMake/3.18.4-GCCcore-10.2.0
module load git/2.28.0-GCCcore-10.2.0-nodocs
module load Python/3.8.6-GCCcore-10.2.0
module load GSL/2.6-GCC-10.2.0
module load Boost/1.74.0-GCC-10.2.0
module load Eigen/3.3.9-GCCcore-10.2.0
module load pybind11/2.6.0-GCCcore-10.2.0
module load Clang/11.0.1-GCCcore-10.2.0
module load PLINK/2.00a2.3_x86_64
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
# Convert files into hacky phased plink1 data
threads files --hap_gz snp_data.hap.gz --sample snp_data.sample --bfile snp_data
# Run the inference
threads infer --bfile snp_data --map_gz snp_data.map.gz --mutation_rate 1.65e-8 --demography CEU.demo --modality snp --threads snp_data.threads
# Convert to .argn and/or .tsz
threads convert --threads snp_data.threads --argn snp_data.argn --tsz snp_data.tsz

# And for sequencing data:
threads files --hap_gz seq_data.hap.gz --sample seq_data.sample --bfile seq_data
threads infer --bfile seq_data --map_gz seq_data.map.gz --mutation_rate 1.65e-8 --demography CEU.demo --modality seq --threads seq_data.threads
threads convert --threads seq_data.threads --argn seq_data.argn --tsz seq_data.tsz
