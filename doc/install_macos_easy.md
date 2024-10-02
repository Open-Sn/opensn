# Easy Install on MacOS

The following instructions were tested on MacOS Sonoma (= version 14).

## Step 1 - Install Development Tools

The following packages should be installed via [Homebrew](https://brew.sh):

- make
- cmake
- gcc
- curl
- bison
- m4
- doxygen
- sphinx

To install any missing packages use
```shell
brew install <package name>
```
Pay close attention to the output of `brew`. It may be necessary to add
locations to your `PATH` environment variable to use the newly installed
utilities.

**Important:** The `sphinx` package may not be available via `Homebrew`. It
is only needed to build the OpenSn documentation.

## Step 2 - Install MPICH

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
If successful, set the following environment variables in the `~/.bashrc` or
`~/.bash_profile` script:
```shell
export MPI_DIR="/path/to/mpich/install"
export PATH="${MPI_DIR}/bin:${PATH}"

export CC="${MPI_DIR}/mpicc"
export CXX="${MPI_DIR}/mpicxx"
export FC="${MPI_DIR}/mpifort"
export F77="${MPI_DIR}/mpif77"
```

## Step 3 - Clone OpenSn

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

## Step 4 - Set Up the Environment

**Important:** XCode 15's linker breaks a number of things in the name of
progress. It may be necessary to modify `configure_dependencies.py`
to use PETSc 3.20.x. It may also require that you add the line
`link_libraries("-ld_classic")` to the OpenSn `CMakeLists.txt`
file.

Next, run the script to compile the necessary dependencies with
```shell
mkdir -p /path/to/dependencies
cd /path/to/opensn
python3 tools/configure_dependencies.py -d /path/to/dependencies
```
It is recommended that `path/to/dependencies` be outside the OpenSn
source tree. Set environment variables for building OpenSn:
```shell
export CMAKE_PREFIX_PATH=/path/to/dependencies:$CMAKE_PREFIX_PATH`
```

## Step 5 - Configure and Build OpenSn

OpenSn is configured within a build directory with
```shell
cd /path/to/opensn
mkdir build
cd build
cmake ..
```
To configure with support for building the documentation use
```shell
cd /path/to/opensn
mkdir build
cd build
cmake -DOPENSN_WITH_DOCS=ON ..
```
In general, the build directory will be within the source tree.
Once configuration is complete, OpenSn can then be built within
the build directory via
```shell
make -j<N>
```

Note: OpenSn may need to be reconfigured with dependency changes, the addition
of new files, etc. When this occurs, clear the `build` directory and repeat
the configuration process above.

## Step 6 - Run Regression Tests

To check if the code compiled correctly execute the test scripts:
```shell
cd /path/to/opensn
test/run_tests -j<N>
```

## Step 7 - OpenSn Documentation

If you configured the OpenSn build environment with support for building the
documentation (see **Step 5**), these instructions will help you install the
necessary tools and build the documentation.

To generate the documentation from your local working copy of OpenSn, you need
to use `pip3` to install the required Python packages:
```bash
pip3 install breathe myst-parser sphinx_rtd_theme
```

Then, from your `build` directory, you can run the command `make doc` to generate
the documentation:
```bash
cd build
make doc
```

Once the build process has completed, you can view the generated documentation by
opening
`opensn/build/doc/index.html` in your favorite web browser.
