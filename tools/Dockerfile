FROM ubuntu:24.04

ARG NUM_CORES=1

RUN apt update -y && \
    apt-get install -y build-essential \
    python3-dev \
    python3-pip \
    python3-venv \
    git \
    cmake \
    libopenmpi-dev \
    flex \
    doxygen \
    pybind11-dev \
    pandoc \
    curl

RUN python3 -m venv /venv && \
    . /venv/bin/activate && \
    git clone https://github.com/Open-Sn/opensn.git

ENV PATH="/venv/bin:$PATH"

RUN cd /opensn/tools && \
    mkdir build_deps && \
    cd build_deps && \
    cmake -DCMAKE_INSTALL_PREFIX=../tpls ../dependencies && \
    make -j ${NUM_CORES} && \
    cd .. && \
    rm -rf build_deps

RUN cd /opensn/doc && \
    pip install -r requirements.txt

ENV CMAKE_PREFIX_PATH=/opensn/tools/tpls

RUN cd /opensn && \
    mkdir -p build && \
    cd build && \
    cmake .. && \
    make -j ${NUM_CORES}

RUN cd /opensn && \
    pip install .

# Common Python packages for data analysis and visualization
RUN pip install matplotlib pandas numpy scipy ipython jupyter

RUN mkdir /work
