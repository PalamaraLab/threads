FROM ubuntu:22.04
RUN apt-get update
RUN apt-get install -y \
    python3 \
    python3-venv \
    python3-pip \
    git \
    clang \
    cmake \
    python3-dev \
    libboost-dev \
    libeigen3-dev \
    pybind11-dev

RUN pip install --upgrade pip setuptools wheel
RUN pip install cmake ninja
