# Compiling on Personal Linux Machines

The following instructions were tested on Ubuntu 18.04 LTS; other Linux
distributions might require some minor tweaking.

Some portions of this document indicate the use of `sudo`. If you do not have 
administrative privileges then have you system administrator assist you in 
these steps.

### Step 1 - Installing GCC, GFortran and the Basic Environment

GCC is used to build and install Chi-Tech.
GFortran and Python is used during the installation of PETSc
(which Chi-Tech uses as a linear algebra backend) and
OpenGL is required by VTK (used by Chi-Tech for visualization).
These packages will therefore also need to be installed.

Check to see if gcc is installed

```bash
gcc --version
```

This should display something like this:

    $ gcc --version
    gcc (Ubuntu 11.3.0-1ubuntu1~22.04.1) 11.3.0
    Copyright (C) 2021 Free Software Foundation, Inc.
    This is free software; see the source for copying conditions.  There is NO
    warranty; not even for MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.

If this is not displayed then install gcc using the following command:

```bash
sudo apt-get install build-essential gfortran
```

If you still do not get the appropriate version when running either ``gcc --version``
or ``gfortran --version`` then there are many online resources which may be able to
assist you.

Now, install the remaining packages needed to build Chi-Tech and its dependencies:

```bash
sudo apt-get install cmake python git zlib1g-dev libx11-dev unzip
```

### Step 2 - An MPI Flavor

Install either OpenMPI or MPICH. If you have MOOSE or deal.ii installed then you
are probably already set to go. **MPICH** is recommended for better performance.

```bash
sudo apt-get install mpich
```

To check if this is working properly, simply do the mpi version of checking for
a C++ compiler:

```bash
mpicc --version
```

Which should display the same message the gcc call did, i.e.

    $ mpicc --version
    gcc (Ubuntu 11.3.0-1ubuntu1~22.04.1) 11.3.0
    Copyright (C) 2021 Free Software Foundation, Inc.
    This is free software; see the source for copying conditions.  There is NO
    warranty; not even for MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.

In order to ensure that all installations downstream of this install with the 
required mpi wrappers, do the following:
```shell
export CC=mpicc
export CXX=mpicxx
export FC=mpifort
```

### Step 3 - PETSc

The current supported version is
[petsc version 3.17.0](https://ftp.mcs.anl.gov/pub/petsc/release-snapshots/petsc-3.17.0.tar.gz).

Return to your *projects* folder (or whatever you chose to place stuff). Do
the following. Download and extract the archive file to a folder of your choice
then navigate into the directory containing the "configure" script.

```bash
wget https://ftp.mcs.anl.gov/pub/petsc/release-snapshots/petsc-3.17.0.tar.gz
tar -zxf petsc-3.17.0.tar.gz
cd petsc-3.17.0
```

and execute the following:

```bash
./configure  \
--prefix=$PWD/../petsc-3.17.0-install  \
--download-hypre=1  \
--with-ssl=0  \
--with-debugging=0  \
--with-pic=1  \
--with-shared-libraries=1  \
--download-fblaslapack=1  \
--download-metis=1  \
--download-parmetis=1  \
--download-superlu_dist=1  \
--with-cxx-dialect=C++11  \
--with-64-bit-indices \
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
that you should use. Execute this command.

After a successful make PETSc should indicate the command used to install
the header files and libraries. Execute the associated command as well.

To test whether the system has been installed correctly execute:

```bash
make test
```

The final step is to generate the environment variable for the install
directory of PETSc:

```bash
export PETSC_ROOT=$PWD/install
```

Again, this is also something you'd like to add to your bash profile.

### Step 4 - Install the Visualization Tool Kit

In your projects folder install VTK using the following commands:

```bash
mkdir VTK
cd VTK
wget https://www.vtk.org/files/release/9.1/VTK-9.1.0.tar.gz
tar -zxf VTK-9.1.0.tar.gz
cd VTK-9.1.0
mkdir build
cd build
cmake -DCMAKE_INSTALL_PREFIX=$PWD/../install  + \
-DBUILD_SHARED_LIBS:BOOL=ON  + \
-DVTK_Group_MPI:BOOL=ON  + \
-DVTK_GROUP_ENABLE_Qt=NO  + \
-DVTK_GROUP_ENABLE_Rendering=NO  + \
-DVTK_GROUP_ENABLE_Imaging=NO  + \
-DVTK_GROUP_ENABLE_StandAlone=WANT  + \
-DVTK_GROUP_ENABLE_Web=NO  + \
-DVTK_BUILD_TESTING:BOOL=OFF  + \
-DCMAKE_BUILD_TYPE=Release  + \
-DCMAKE_CXX_FLAGS=-std=c++11  + \
 ../
```

This will take a while after which you need to execute:

```bash
make -j4 && make install
```

This will also take a while. Finally, export the *VTK_DIR* variable:

```bash
export VTK_DIR=$PWD/install
```

Again, this is also something you'd like to add to your bash profile.

### Step 5 - Install Lua

#### Step 5a - Install ncurses
Navigate to a folder to install `ncurses`
```shell
wget https://invisible-mirror.net/archives/ncurses/ncurses-6.1.tar.gz
tar -zxf ncurses-6.1.tar.gz
cd ncurses-6.1
./configure --prefix=$PWD/chi_build
```

After the configure completes then
```shell
make -j4
```
and then
```shell
make install
```

#### Step 5b - Install readline
Navigate to a folder to install `readline`
```shell
wget ftp://ftp.gnu.org/gnu/readline/readline-8.0.tar.gz
tar -zxf readline-8.0.tar.gz
cd readline-8.0
./configure --prefix=$PWD/chi_build
```

After the configure completes then
```shell
make -j4
```
and then
```shell
make install
```

#### Step 5c - The actual lua installation

Download and extract **Lua** version 5.3.5 from https://www.lua.org. 
Before installing **Lua** edit the Makefile and set INSTALL_TOP to your desired
install location.  Install **Lua** as follows:
```bash
$ export CPATH=<path-to-readline-include>:$CPATH
$ export LIBRARY_PATH=<path-to-readline-lib>:<path-to-ncurses-lib>:$LIBRARY_PATH
$ make linux
$ make local
```
If the install complains about missing **readline** includes or libraries, it may
be necessary to install **readline** first, that would probably require 
**ncurses** in turn.

Set the LUA_ROOT environment variable to the **Lua** install location:
```bash
    $export LUA_ROOT=/Path/to/Lua
```

Add the export command to your bash profile.

### Step 6 - Build OpenSn

Clone the **OpenSn** repository.  Go the folder where you want to keep OpenSn relevant stuff:
```bash
    $ git clone https://github.com/Open-Sn/opensn.git /path/to/opensn
```

Go to the opensn folder and type:
```bash
    $ cd opensn
    $ ./configure.sh
```
The configure script will generate the CMake build scripts.

In the main build directory (e.g., `open/build/`), execute:
```bash
    $ make -j4
```

You can also use -j8 even if you don't have 8 processors, the make command
will use threading where possible.

### Step 7 - Run Regression Tests

To check if the code compiled correctly, execute the test scripts:

```bash
    $ test/run_tests -d test/ -j8
```


### Step 8 - OpenSn Documentation

You can either access the documentation online [here](https://xxx.io), or generate it locally.

To generate the documentation from your local working copy, first make sure
Doxygen and LaTeX are installed:

```bash
sudo apt-get install doxygen texlive
```

The documentation is contained in the *doc* folder and can be generated
using a script provided in that folder:

```bash
cd doc/
./YReGenerateDocumentation.sh
```

Once finished, you can view the generated documentation by opening

```bash
doc/HTMLdocs/html/index.html
```

in your favorite browser.