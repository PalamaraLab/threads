FROM ubuntu:22.04

# Base python dependencies
RUN apt-get update
RUN apt-get install -y \
    python3 \
    python3-venv \
    python3-pip

# Packages for local testing
RUN pip install --upgrade pip setuptools wheel
RUN pip install pytest h5py numpy

# Packages for local build
RUN pip install cmake ninja

# Dev tools require an additional update
RUN apt-get update
RUN apt-get install -y \
    git \
    gdb \
    lcov

# Dev dependencies
RUN apt-get install -y \
    python3-dev \
    libboost-dev \
    pybind11-dev
