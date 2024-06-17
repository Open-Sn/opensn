# Easy Install on Linux Machines

The following instructions were tested on Ubuntu 22.04.4 LTS. Other Linux
distributions might require some minor tweaking.

Some portions of these instructions require the use of `sudo`. If you do not have
administrative privileges, have your system administrator assist you with these
steps.

Note that these instructions assume you are using the `bash` shell and may differ
slightly for other shells.

## Step 1 - Install Development Tools

The following packages are required for OpenSn installation and development:

1. A recent version of clang++/g++ that supports C++17
2. gfortran (required by the BLAS component of PETSc)
3. flex and bison (required by the PTSCOTCH component of PETSc)
4. Python 3 v3.9+ and pip (required by PETSc and OpenSn)
5. ncurses v5 (required by Lua)
6. Git version control system
7. CMake v3.12+
8. MPI (OpenMPI, MPICH, and MVAPICH have been tested)
9. Doxygen and Sphinx (required for generating the OpenSn documentation)
10. curl (required to download and install the OpenSn third-party dependencies)

Most of these packages can be installed using the package manager available with
your Linux distribution (e.g., `apt`, `yum`, etc.). For example, on Ubuntu, you
can use the following command to install all of these packages:

```bash
$ sudo apt install build-essential gfortran python3 \
git cmake libopenmpi-dev flex bison \
libncurses5-dev python3-pip doxygen sphinx curl
```

To ensure that all third-party packages use the required MPI compiler wrappers,
the following environment variables must be set:

```bash
    $ export CC=mpicc
    $ export CXX=mpicxx
    $ export FC=mpifort
```

## Step 2 - Clone OpenSn

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

## Step 3 - Install Third-Party Libraries

**Important:** We recommend creating a separate directory for building the
**OpenSn** dependencies.

Assuming you created a directory named `dependencies` to be used for building the
required third-party packages, the following command automates the install and
build process:

```bash
    $ cd opensn
    $ python3 resources/configure_dependencies.py -d ../path/to/dependencies/directory
```

## Step 4 - Configure Environment

Before compiling **OpenSn**, you must add the location of the third-party
libraries to your `CMAKE_PREFIX_PATH` environment variable. We have provided a
script to do this for you. You can update your environment variables by running
the command:

```bash
    $ source ../dependencies/configure_deproots.sh
```
**Important:** It may be a good idea to source this script int your `.bashrc`
file so that you don't need to specify the path every time you need to re-run
`cmake`.

## Step 5 - Build OpenSn

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

For more information on building the documentation, see **Step 7** below.

## Step 6 - Run Regression Tests

To run the regression tests, simply run `make test` from the build directory.
This will run all of the regression tests in the `opensn/test` directory.

## Step 7 - Build the OpenSn Documentation

If you configured the **OpenSn** build environment with support for building the
documentation (see **Step 5**), these instructions will help you install the
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