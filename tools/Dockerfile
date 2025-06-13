FROM ubuntu:24.04

ARG NUM_CORES=1

ENV CC=/usr/bin/mpicc
ENV CXX=/usr/bin/mpicxx
ENV FC=/usr/bin/mpifort

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
    python3 configure_dependencies.py -d dependencies -j ${NUM_CORES}

ENV CMAKE_PREFIX_PATH=/opensn/tools/dependencies

RUN cd /opensn && \
    mkdir -p build && \
    cd build && \
    cmake .. && \
    make -j ${NUM_CORES}

RUN cd /opensn && \
    pip install .

# Common Python packages for data analysis and visualization
RUN pip install matplotlib pandas numpy ipython jupyter

RUN mkdir /work
