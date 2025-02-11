# Install on MacOS

The following instructions assume that all the prerequisite software packages
will be installed by hand.  These packages may also be installed via Homebrew
or MacPorts.

## Step 1 - Install GCC
```shell
brew install gcc
```

This will generally install most recent version of `gcc` to `/usr/local/bin`.
In this directory, there will be links such as `gcc-13`, `g++-13`, and
`gfortran-13` pointing to the Homebrew installation.

Note: `gfortran` will be needed by PETSc to build its dependencies.

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

## Step 3 - Install PETSc

Download a suitable version of
[PETSc](https://web.cels.anl.gov/projects/petsc/download/release-snapshots/).
The current supported version is
[PETSc 3.17.0](https://web.cels.anl.gov/projects/petsc/download/release-snapshots/petsc-3.17.0.tar.gz).
Unpack the source code into `/path/to/petsc`.
Configure PETSc with the appropriate dependencies:
```shell
cd /path/to/petsc
./configure  \
--prefix=/path/to/petsc/install \
--with-shared-libraries=1  \
--with-ssl=0  \
--with-debugging=0  \
--with-pic=1  \
--with-64-bit-indices=1 \
--download-hypre=1  \
--download-fblaslapack=1  \
--download-metis=1  \
--download-parmetis=1  \
--download-superlu_dist=1  \
--download-ptscotch=1  \
CC=$CC CXX=$CXX FC=$FC \
CFLAGS='-fPIC -fopenmp'  \
CXXFLAGS='-fPIC -fopenmp'  \
FFLAGS='-fPIC -fopenmp'  \
FCFLAGS='-fPIC -fopenmp'  \
F90FLAGS='-fPIC -fopenmp'  \
F77FLAGS='-fPIC -fopenmp'  \
COPTFLAGS='-O3 -march=native -mtune=native'  \
CXXOPTFLAGS='-O3 -march=native -mtune=native'  \
FOPTFLAGS='-O3 -march=native -mtune=native'  \
PETSC_DIR=$PWD
```
If the configuration fails then consult PETSc's user documentation.

Upon completion of the configure step, PETSc will provide the make command
that you should use. The addition of the build option `OMAKE_PRINTDIR=make`
may be required in some circumstances with GNU compilers.

After a successful build, PETSc will provide instructions for installing and
checking the installations. Follow these instructions.

After a successful install, add the following environment variables to the
`~/.bashrc` or `~/.bash_profile` script:
```shell
export PETSC_ROOT="/path/to/petsc/install"
export CMAKE_PREFIX_PATH="${PETSC_ROOT}${CMAKE_PREFIX_PATH:+:${CMAKE_PREFIX_PATH}}"
```

## Step 4 - Install the Visualization Tool Kit

Download [VTK 9.1.0]( https://www.vtk.org/files/release/9.1/VTK-9.1.0.tar.gz) or
[VTK 9.3.0](https://www.vtk.org/files/release/9.3/VTK-9.3.0.tar.gz)
into a suitable location and upack it into `/path/to/vtk`. Assuming VTK will
be built in `/path/to/vtk/build` and installed in `/path/to/vtk/install`, execute
the following to configure, build, and install VTK:
```shell
mkdir -p /path/to/vtk/build
cd /path/to/vtk/build

cmake \
-DCMAKE_INSTALL_PREFIX=/path/to/vtk/install \
-DBUILD_SHARED_LIBS=ON \
-DVTK_USE_MPI=ON \
-DVTK_GROUP_ENABLE_StandAlone=WANT \
-DVTK_GROUP_ENABLE_Rendering=DONT_WANT \
-DVTK_GROUP_ENABLE_Imaging=DONT_WANT \
-DVTK_GROUP_ENABLE_Web=DONT_WANT \
-DVTK_GROUP_ENABLE_Qt=DONT_WANT \
-DVTK_MODULE_USE_EXTERNAL_VTK_hdf5=ON \
-DCMAKE_BUILD_TYPE=Release \
/path/to/vtk

make -j<N>
make install
```

After a successful install, set the following environment variables in the
`~/.bashrc` or `~/.bash_profile` script:
```shell
export VTK_DIR="/path/to/vtk/install"
export CMAKE_PREFIX_PATH="${VTK_DIR}:${CMAKE_PREFIX_PATH}"
```

## Step 5 - Install Lua

Lua requires the
[readline](ftp://ftp.gnu.org/gnu/readline/readline-8.0.tar.gz) and
[ncurses](https://invisible-mirror.net/archives/ncurses/ncurses-6.1.tar.gz)
packages. If these are not installed, either install them from source or
using Homebrew with:
```shell
brew install readline ncurses
```
Once installed, set the following environment variables in the `~/.bashrc`
or `~/.bash_profile` script:
```shell
export LIBRARY_PATH="/path/to/readline/install/lib:${LIBRARY_PATH}"
export LIBRARY_PATH="/path/to/ncurses/install/lib:${LIBRARY_PATH}"
export CPATH="/path/to/readline/install/include:${CPATH}"
export CPATH="/path/to/ncurses/install/include:${CPATH}"
```

Download [Lua 5.4+](https://www.lua.org/ftp/) into a suitable location and
unpack it into `path/to/lua`.
[Lua 5.4.6](https://www.lua.org/ftp/lua-5.4.6.tar.gz) is recommended.
Execute the following to build and install Lua:
```shell
cd /path/to/lua
make macosx CC=gcc-<version> MYCFLAGS=-fPIC MYLIBS=-lncurses -j<N>
make install INSTALL_TOP=/path/to/lua/install
```
where `<version>` is the GCC version installed via Homebrew.

After a successful installation, set the following environment variables in
the `~/.bashrc` or `~/.bash_profile` script:
```shell
export LUA_ROOT="/path/to/lua/install"
export CMAKE_PREFIX_PATH="${LUA_ROOT}:${CMAKE_PREFIX_PATH}"
```

## Step 6 - Install Caliper

Download [Caliper 2.10.0](https://github.com/LLNL/Caliper/archive/refs/tags/v2.10.0.tar.gz) or
into a suitable location and upack it into `/path/to/caliper`. Assuming Caliper will
be built in `/path/to/caliper/build` and installed in `/path/to/caliper/install`, execute
the following to configure, build, and install Caliper:
```shell
mkdir -p /path/to/caliper/build
cd /path/to/caliper/build

cmake \
-DCMAKE_INSTALL_PREFIX=/path/to/caliper/install \
-DWITH_MPI=ON \
-DWITH_KOKKOS=OFF \
/path/to/caliper

make -j<N>
make install

After a successful install, set the following environment variables in the
`~/.bashrc` or `~/.bash_profile` script:
```shell
export CALIPER_DIR="/path/to/caliper/install"
export CMAKE_PREFIX_PATH="${CALIPER_DIR}:${CMAKE_PREFIX_PATH}"
```

## Step 7 - Clone OpenSn

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

## Step 8 - Configure and Build OpenSn

OpenSn is configured within a build directory with
```shell
cd /path/to/opensn
mkdir build
cd build
cmake ..
```
This will configure the project for building it.
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

## Step 9 - Run Regression Tests

To check if the code compiled correctly execute the test scripts:
```shell
cd /path/to/opensn
test/run_tests -j<N>
build/test/opensn-unit
```

## Step 10 - OpenSn Documentation

If you configured the OpenSn build environment with support for building the
documentation (see **Step 8**), these instructions will help you install the
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
