FROM ubuntu:24.04

ENV CC=/usr/bin/mpicc
ENV CXX=/usr/bin/mpicxx
ENV FC=/usr/bin/mpifort

ENV NPROC=20

RUN apt update -y && \
    apt-get install -y build-essential \
    gfortran \
    python3-dev \
    python3-pip \
    python3-venv \
    git \
    cmake \
    libopenmpi-dev \
    flex \
    bison \
    libhdf5-mpi-dev \
    doxygen \
    pybind11-dev \
    curl

RUN python3 -m venv /venv && \
    . /venv/bin/activate && \
    git clone https://github.com/Open-Sn/opensn.git

ENV PATH="/venv/bin:$PATH"

RUN cd /opensn/tools && \
    mkdir dependencies && \
    python3 configure_dependencies.py -d dependencies

ENV CMAKE_PREFIX_PATH=/opensn/tools/dependencies

RUN cd /opensn && \
    mkdir -p build && \
    cd build && \
    cmake .. && \
    make -j20

RUN cd /opensn && \
    pip install .

RUN pip install matplotlib pandas numpy ipython

RUN mkdir /work
