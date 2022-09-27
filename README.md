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
cd ..
```