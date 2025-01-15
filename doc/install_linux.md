# Install on Linux Machines

The following instructions were tested on Ubuntu 22.04.4 LTS. Other Linux
distributions might require some minor tweaking.

Some portions of these instructions require the use of `sudo`. If you do not
have administrative privileges, have your system administrator assist you with
these steps.

Note that these instructions assume you are using the `bash` shell and may
differ slightly for other shells.

## Step 1 - Install Development Tools

The following packages are required for **OpenSn** installation and development:

1. A recent version of clang++/g++ that supports C++17
2. gfortran (required by the BLAS component of PETSc)
3. flex and bison (required by the PTSCOTCH component of PETSc)
4. Python 3 v3.9+ and pip (required by PETSc and OpenSn)
5. ncurses v5 (required by Lua)
6. Git version control system
7. CMake v3.12+
8. MPI (OpenMPI, MPICH, and MVAPICH have been tested)
9. Doxygen and Sphinx (required for generating the OpenSn documentation)

Most of these packages can be installed using the package manager available with
your Linux distribution (e.g., `apt`, `yum`, etc.). For example, on Ubuntu, you
can use the following command to install all of these packages:

```bash
$ sudo apt install build-essential gfortran python3 \
git cmake libopenmpi-dev flex bison \
libncurses5-dev python3-pip doxygen sphinx
```

To ensure that all third-party packages use the required MPI compiler wrappers,
the following environment variables must be set:

```bash
    $ export CC=mpicc
    $ export CXX=mpicxx
    $ export FC=mpifort
```

**Important:** We recommend creating a separate directory for building the
**OpenSn** dependencies. 

The following steps assume that you have created a directory named `dependencies`
to be used for compiling and installing the required third-party libraries.

## Step 2 - Install PETSc

The current supported version is
[**PETSc** version 3.17.0+](https://ftp.mcs.anl.gov/pub/petsc/release-snapshots/petsc-3.17.0.tar.gz).

In your `dependencies` directory, install **PETSc** using the following commands:

```bash
wget https://ftp.mcs.anl.gov/pub/petsc/release-snapshots/petsc-3.17.0.tar.gz
tar -zxf petsc-3.17.0.tar.gz
cd petsc-3.17.0
```

To configure **PETSc** with the required build options, run:

```bash
./configure  \
--prefix=/path/to/dependencies/directory  \
--download-hypre=1  \
--with-ssl=0  \
--with-debugging=0  \
--with-pic=1  \
--with-shared-libraries=1  \
--download-fblaslapack=1  \
--download-metis=1  \
--download-parmetis=1  \
--download-superlu_dist=1  \
--download-ptscotch=1  \
--with-cxx-dialect=C++11  \
--with-64-bit-indices \
CC=$CC CXX=$CXX FC=$FC \
CFLAGS='-fPIC -fopenmp'  \
CXXFLAGS='-fPIC -fopenmp'  \
FFLAGS='-fPIC -fopenmp'  \
FCFLAGS='-fPIC -fopenmp'  \
F90FLAGS='-fPIC -fopenmp'  \
F77FLAGS='-fPIC -fopenmp'  \
COPTFLAGS='-O3'  \
CXXOPTFLAGS='-O3'  \
FOPTFLAGS='-O3'  \
PETSC_DIR=$PWD
```

If the configuration fails, consult **PETSc**'s user documentation. 

Follow the **PETSc** build prompts to complete the **PETSc** installation.

## Step 3 - Install the Visualization Tool Kit (VTK)

In your `dependencies` directory, install **VTK** using the following commands:

```bash
mkdir VTK
cd VTK
wget https://www.vtk.org/files/release/9.3/VTK-9.3.0.tar.gz
tar -zxf VTK-9.3.0.tar.gz
cd VTK-9.3.0
mkdir build
cd build
cmake -DCMAKE_INSTALL_PREFIX=/path/to/dependencies/directory \
-DBUILD_SHARED_LIBS:BOOL=ON \
-DVTK_Group_MPI:BOOL=ON \
-DVTK_GROUP_ENABLE_Qt=NO \
-DVTK_GROUP_ENABLE_Rendering=NO \
-DVTK_GROUP_ENABLE_Imaging=NO \
-DVTK_GROUP_ENABLE_StandAlone=WANT \
-DVTK_GROUP_ENABLE_Web=NO \
-DVTK_MODULE_USE_EXTERNAL_VTK_hdf5=ON \
-DVTK_BUILD_TESTING:BOOL=OFF \
-DCMAKE_BUILD_TYPE=Release \
-DCMAKE_CXX_FLAGS=-std=c++11 \
../
```

After `cmake` has completed configuring the build, run:

```bash
make -j && make install
```

to install **VTK**.

## Step 4 - Install Lua

Download and extract **Lua** version 5.3.6+ from https://www.lua.org to
`dependencies`. Install **Lua** as follows:

```bash
$ make linux MYCFLAGS=-fPIC -j
$ make install INSTALL_TOP=/path/to/dependencies/directory
```

If the install complains about missing **readline** includes or libraries, it
may be necessary to install **readline** and **ncurses**.

## Step 5 - Install Caliper

In your `dependencies` directory, install **Caliper** using the following commands:

```bash
mkdir caliper
cd caliper
wget https://github.com/LLNL/Caliper/archive/refs/tags/v2.10.0.tar.gz
tar -zxf v2.10.0.tar.gz
cd Caliper-2.10.0
mkdir build
cd build
cmake -DCMAKE_INSTALL_PREFIX=/path/to/dependencies/directory  + \
-DWITH_MPI=ON \
-DWITH_KOKKOS=OFF \
 ../
```

After `cmake` has completed configuring the build, run:

```bash
make -j && make install
```

to install **Caliper**.


## Step 6 - Configure Environment

Before compiling **OpenSn**, you must add the location of the third-party
libraries to your `CMAKE_PREFIX_PATH` environment variable. This can be
accomplished with the following command:

```bash
    $ export CMAKE_PREFIX_PATH=/path/to/dependencies${CMAKE_PREFIX_PATH:+:${CMAKE_PREFIX_PATH}}`
```

**Important:** It may be a good idea to add the `CMAKE_PREFIX_PATH` variable to
your `.bashrc` file so that you don't need to specify the path every time you
need to re-run `cmake`.

## Step 7 - Clone OpenSn

**Important:**  If you want to contribute to **OpenSn**, it is strongly
recommended that you first fork the **OpenSn** repository then clone your fork.

To clone the **OpenSn** repository:

```bash
    $ git clone https://github.com/Open-Sn/opensn.git
```

To clone your fork of **OpenSn**:

```bash
    $ git clone https://github.com/<username>/opensn.git
```

## Step 8 - Build OpenSn

To build **OpenSn**, create a build directory in the top-level **OpenSn**
directory and run `cmake` to generate the build files and `make` to compile
**OpenSn**:

```bash
    $ mkdir build
    $ cd build
    $ cmake ..
    $ make -j<N>
```

To configure **OpenSn** for building the documentation, in addition to the
**OpenSn** application, add the `-DOPENSN_WITH_DOCS` option to `cmake`:

```bash
    $ mkdir build
    $ cd build
    $ cmake -DOPENSN_WITH_DOCS=ON ..
    $ make -j<N>
```

For more information on building the documentation, see **Step 10** below.

## Step 9 - Run Regression Tests

To run the regression tests, simply run `make test` from the build directory.
This will run all of the regression tests in the `opensn/test` directory.

## Step 10 - Build the OpenSn Documentation

If you configured the **OpenSn** build environment with support for building the
documentation (see **Step 8**), these instructions will help you install the
necessary tools and build the documentation.

To generate the documentation from your local working copy of **OpenSn**, you
need to use `pip` to install the required **Python** packages:

```bash
pip install breathe myst-parser sphinx_rtd_theme
```

Then, from your `build` directory, you can run the command `make doc` to generate
the documentation:

```bash
cd build
make doc
```

Once the build process has completed, you can view the generated documentation by
opening `opensn/build/doc/index.html` in your favorite web browser.
