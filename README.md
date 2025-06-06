# Threads

This software implements the Threads algorithm, described in

`Á. F. Gunnarsson, J. Zhu, B. C. Zhang, Z. Tsangalidou, A. Allmont, P. Palamara. A scalable approach for genome-wide inference of ancestral recombination graphs. bioRxiv, 2024.`

The user manual for threads can be found [here](https://palamaralab.github.io/software/threads/).

## Rescomp user installation instructions

Clone the threads-arg repo:

```sh
git clone https://github.com/PalamaraLab/threads.git
cd threads
```

Load development modules:

```sh
module load Python/3.11.3-GCCcore-12.3.0
module load Boost/1.82.0-GCC-12.3.0
module load HDF5/1.14.0-gompi-2023a
```

Create a new venv, activate and `pip install` to build:

```sh
python -m venv venv
. venv/bin/activate
pip install --upgrade pip setuptools wheel
pip install .
```

For active development use `-e` and `[dev]` to for additional dependencies:

```sh
pip install -e .[dev]
```

## For developers: making a release

- Bump the version number in [pyproject.toml](pyproject.toml), [CMakeLists.txt](CMakeLists.txt)
- Update [RELEASE_NOTES.md](RELEASE_NOTES.md)
- Push changes and check that all [GitHub workflows](https://github.com/PalamaraLab/threads/actions) pass
- Tag the commit in Git using syntax `vX.Y.Z`
- Make a release on GitHub, which should trigger a new build that will upload Python wheels to PyPI
- If building wheels for a new Python version:
  - Update classifiers in [pyproject.toml](pyproject.toml)
  - Check whether version of `pypa/cibuildwheel` needs updating in [build-wheels.yml](.github/workflows/build-wheels.yml)
