# Easy Install on MacOS

The following instructions were tested on MacOS Monterey (= version 12).

## Prerequisites

The following packages should be installed via Homebrew
- make
- cmake
- gcc
- wget

To install any missing packages use
```shell
brew install <package name>
```

## Step 1 - MPICH

Download a suitable version of [MPICH 4+](https://www.mpich.org/static/downloads).
Versions above 4.0 are recommended.
Unpack the source into `/path/to/mpich`. Assuming MPICH will be built in
`/path/to/mpich/build` and installed in `/path/to/mpich/install`, execute the
following to configure, build, and install MPICH:
```shell
mkdir -p /path/to/mpich/build
cd /path/to/mpich/build

/path/to/mpich/configure \
--prefix=/path/to/mpich/install \
CC=gcc-13 CXX=g++-13 FC=gfortran-13 F77=gfortran-13 \
FCFLAGS=-fallow-argument-mismatch FFLAGS=-fallow-argument-mismatch

make -j<N>
make install
```
It is recommended by MPICH that the build and install directories be outside
the source tree `/path/to/mpich`.

Check the installation:
```shell
/path/to/mpich/install/bin/mpicc --version
```
If the installation was successful, a message similar to
```
gcc-13 (Homebrew GCC 13.2.0) 13.2.0
Copyright (C) 2023 Free Software Foundation, Inc.
This is free software; see the source for copying conditions.  There is NO
warranty; not even for MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
```
should appear. Ensure that the message displays the desired compilers.
If successful, set the following evironment variables in the `~/.bashrc` or
`~/.bash_profile` script:
```shell
export MPI_DIR="/path/to/mpich/install"
export PATH="${MPI_DIR}/bin:${PATH}"

export CC="${MPI_DIR}/mpicc"
export CXX="${MPI_DIR}/mpicxx"
export FC="${MPI_DIR}/mpifort"
export F77="${MPI_DIR}/mpif77"
```

## Step 2 - Clone OpenSn

The remainder of the instructions are identical to
[Easy Linux instructions](./install_ubuntu_easy.md).

Note:  If you want to contribute to OpenSn, it is strongly recommended
to first fork the OpenSn repository into your own Git account and then to
clone your fork.

Clone the OpenSn repository or a fork:
```shell
git clone https://github.com/Open-Sn/opensn.git /path/to/opensn
````
or
```shell
git clone https://github.com/<username>/opensn.git /path/to/opensn
```

## Step 3 - Set Up the Enviroment

Next, run the script to compile the necessary dependencies with
```shell
mkdir -p /path/to/dependencies
cd /path/to/opensn
python3 resources/configure_dependencies.py -d /path/to/dependencies
```
It is recommended that `path/to/dependencies` be outside the OpenSn
source tree. Set environment variables for building OpenSn with the
generated script:
```shell
source /path/to/dependencies/configure_deproots.sh
```

## Step 4 - Configure and Build OpenSn

OpenSn is configured within a build directory with
```shell
cd /path/to/opensn
mkdir build
cd build
cmake ..
```
This will configure the project for building it.
In general, the build directory will be within the source tree.
OpenSn can then be built within the build directory via
```shell
make -j<N>
```

Note: OpenSn may need to be reconfigured with dependency changes, the addition
of new files, etc. When this occurs, clear the `build` directory and repeat
the configuration process above.

## Step 5 - Run Regression Tests

To check if the code compiled correctly execute the test scripts:
```shell
cd /path/to/opensn
test/run_tests -j<N>
```

## Step 6 - OpenSn Documentation

The documentation can be found [online](https://chi-tech.github.io), or
generated locally. To generate the documentation locally, first make sure
doxygen and LaTeX are installed:
```shell
sudo apt-get install doxygen texlive
```

The documentation is contained in the `doc` directory of the OpenSn source
tree and can be generated with
```shell
./YReGenerateDocumentation.sh
```
from within the `doc` directory. Once finished, the generated documentation
can be viewed with
```shell
doc/HTMLdocs/html/index.html
```
in a web browser.