Install on macOS
================

The following instructions assume that all prerequisite software packages will
be installed manually. These packages may also be installed via Homebrew or
MacPorts.

Step 1 - Install GCC
--------------------

To install GCC using Homebrew, run:

.. code-block:: shell

   brew install gcc

This will typically install the latest version of ``gcc`` to ``/usr/local/bin``.
In this directory, you should find symbolic links such as ``gcc-13``,
``g++-13``, and ``gfortran-13`` pointing to the Homebrew installation.

.. note::

   ``gfortran`` is required by PETSc to build its dependencies.

Step 2 - Install MPICH
----------------------

Download a suitable version of `MPICH 4+ <https://www.mpich.org/static/downloads>`_.
Versions above 4.0 are recommended.

Unpack the source into ``/path/to/mpich``. Assuming MPICH will be built in
``/path/to/mpich/build`` and installed in ``/path/to/mpich/install``, execute
the following commands to configure, build, and install MPICH:

.. code-block:: shell

   mkdir -p /path/to/mpich/build
   cd /path/to/mpich/build

   /path/to/mpich/configure \
      --prefix=/path/to/mpich/install \
      CC=gcc-13 CXX=g++-13 FC=gfortran-13 F77=gfortran-13 \
      FCFLAGS=-fallow-argument-mismatch FFLAGS=-fallow-argument-mismatch

   make -j<N>
   make install

.. note::

   MPICH recommends placing the build and install directories outside the source
   tree (``/path/to/mpich``).

To verify the installation, run:

.. code-block:: shell

   /path/to/mpich/install/bin/mpicc --version

If the installation was successful, you should see output similar to:

.. code-block:: text

   gcc-13 (Homebrew GCC 13.2.0) 13.2.0
   Copyright (C) 2023 Free Software Foundation, Inc.
   This is free software; see the source for copying conditions.  There is NO
   warranty; not even for MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.

Ensure the output reflects the correct compilers.

If everything is working correctly, add the following environment variables to
your ``~/.bashrc`` or ``~/.bash_profile`` script:

.. code-block:: shell

   export MPI_DIR="/path/to/mpich/install"
   export PATH="${MPI_DIR}/bin:${PATH}"

   export CC="${MPI_DIR}/mpicc"
   export CXX="${MPI_DIR}/mpicxx"
   export FC="${MPI_DIR}/mpifort"
   export F77="${MPI_DIR}/mpif77"

Step 3 - Install PETSc
-----------------------

Download a suitable version of `PETSc <https://web.cels.anl.gov/projects/petsc/download/release-snapshots/>`_.  
The currently supported version is `PETSc 3.17.0 <https://web.cels.anl.gov/projects/petsc/download/release-snapshots/petsc-3.17.0.tar.gz>`_.

Unpack the source code into ``/path/to/petsc`` and configure PETSc with the
appropriate dependencies:

.. code-block:: shell

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

If the configuration fails, consult PETSc's user documentation.

After the configure step, PETSc will provide a ``make`` command for building.
You may need to append the build option ``OMAKE_PRINTDIR=make`` in some cases
when using GNU compilers.

Once built successfully, PETSc will also provide instructions for installing and
verifying the installation. Follow those instructions.

After a successful install, add the following environment variables to your
``~/.bashrc`` or ``~/.bash_profile`` script:

.. code-block:: shell

   export PETSC_ROOT="/path/to/petsc/install"
   export CMAKE_PREFIX_PATH="${PETSC_ROOT}${CMAKE_PREFIX_PATH:+:${CMAKE_PREFIX_PATH}}"

Step 4 - Install the Visualization Tool Kit
-------------------------------------------

Download either `VTK 9.1.0 <https://www.vtk.org/files/release/9.1/VTK-9.1.0.tar.gz>`_
or `VTK 9.3.0 <https://www.vtk.org/files/release/9.3/VTK-9.3.0.tar.gz>`_ into a
suitable location and unpack it into ``/path/to/vtk``.

Assuming VTK will be built in ``/path/to/vtk/build`` and installed in
``/path/to/vtk/install``, run the following to configure, build, and install:

.. code-block:: shell

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

After a successful install, set the following environment variables in your
``~/.bashrc`` or ``~/.bash_profile`` script:

.. code-block:: shell

   export VTK_DIR="/path/to/vtk/install"
   export CMAKE_PREFIX_PATH="${VTK_DIR}:${CMAKE_PREFIX_PATH}"

Step 5 - Install Lua
--------------------

Lua requires the `readline <ftp://ftp.gnu.org/gnu/readline/readline-8.0.tar.gz>`_
and `ncurses <https://invisible-mirror.net/archives/ncurses/ncurses-6.1.tar.gz>`_
packages. If these are not installed, you can either install them from source or
use Homebrew:

.. code-block:: shell

   brew install readline ncurses

Once installed, add the following environment variables to your ``~/.bashrc`` or
``~/.bash_profile`` script:

.. code-block:: shell

   export LIBRARY_PATH="/path/to/readline/install/lib:${LIBRARY_PATH}"
   export LIBRARY_PATH="/path/to/ncurses/install/lib:${LIBRARY_PATH}"
   export CPATH="/path/to/readline/install/include:${CPATH}"
   export CPATH="/path/to/ncurses/install/include:${CPATH}"

Download `Lua 5.4+ <https://www.lua.org/ftp/>`_ into a suitable location and
unpack it to ``/path/to/lua``. We recommend using `Lua 5.4.6 <https://www.lua.org/ftp/lua-5.4.6.tar.gz>`_.

To build and install Lua:

.. code-block:: shell

   cd /path/to/lua
   make macosx CC=gcc-<version> MYCFLAGS=-fPIC MYLIBS=-lncurses -j<N>
   make install INSTALL_TOP=/path/to/lua/install

Replace ``<version>`` with the GCC version you installed via Homebrew.

After a successful installation, set the following environment variables in your
``~/.bashrc`` or ``~/.bash_profile`` script:

.. code-block:: shell

   export LUA_ROOT="/path/to/lua/install"
   export CMAKE_PREFIX_PATH="${LUA_ROOT}:${CMAKE_PREFIX_PATH}"

Step 6 - Install Caliper
------------------------

Download `Caliper 2.10.0 <https://github.com/LLNL/Caliper/archive/refs/tags/v2.10.0.tar.gz>`_
into a suitable location and unpack it into ``/path/to/caliper``.

Assuming Caliper will be built in ``/path/to/caliper/build`` and installed in
``/path/to/caliper/install``, run the following:

.. code-block:: shell

   mkdir -p /path/to/caliper/build
   cd /path/to/caliper/build

   cmake \
      -DCMAKE_INSTALL_PREFIX=/path/to/caliper/install \
      -WITH_MPI=ON \
      -WITH_KOKKOS=OFF \
      /path/to/caliper

   make -j<N>
   make install

After a successful install, set the following environment variables in your
``~/.bashrc`` or ``~/.bash_profile`` script:

.. code-block:: shell

   export CALIPER_DIR="/path/to/caliper/install"
   export CMAKE_PREFIX_PATH="${CALIPER_DIR}:${CMAKE_PREFIX_PATH}"

Step 7 - Clone OpenSn
---------------------

.. note::

   If you want to contribute to **OpenSn**, it is strongly recommended that you
   first fork the OpenSn repository into your own Git account and then clone
   your fork.

To clone the OpenSn repository:

.. code-block:: shell

   git clone https://github.com/Open-Sn/opensn.git /path/to/opensn

Or, to clone your fork of OpenSn:

.. code-block:: shell

   git clone https://github.com/<username>/opensn.git /path/to/opensn

Step 8 - Configure and Build OpenSn
-----------------------------------

Lua interface
^^^^^^^^^^^^^

OpenSn is configured within a build directory:

.. code-block:: shell

   cd /path/to/opensn
   mkdir build
   cd build
   cmake ..

This configures the project for building.

To configure OpenSn with support for building the documentation, use:

.. code-block:: shell

   cd /path/to/opensn
   mkdir build
   cd build
   cmake -DOPENSN_WITH_DOCS=ON ..

In general, the build directory is located within the source tree.

Once configuration is complete, OpenSn can be built within the build directory
using:

.. code-block:: shell

   make -j<N>

.. note::

   OpenSn may need to be reconfigured when dependencies change, new files are
   added, etc. In such cases, clear the ``build`` directory and repeat the
   configuration process.

Python console/interface
^^^^^^^^^^^^^^^^^^^^^^^^

**OpenSn** also provides a Python interface. It is available in two formats: a
console application ``opensn`` and a Python module ``pyopensn``.

Classes and functions in the Python interface are detailed in :ref:`pyapi`.

.. attention::

   The console and the module are **not compatible** with each other. Attempting
   to import the module within the console will result in an import error. Users
   should select one approach and maintain consistent coding style throughout.

To compile the console application:

.. code-block:: bash

   mkdir build
   cd build
   cmake -DOPENSN_WITH_PYTHON=ON ..
   make -j<N>

.. danger::

   In the console application, all classes and functions are implicitly imported
   into the ``__main__`` module at startup. Therefore, omit submodule prefixes
   when referring to class or function names. Additionally, avoid redefining any
   **OpenSn** class or function names to prevent naming conflicts.

To compile the module and install in the Python ``site-packages`` path:

.. code-block:: bash

   pip3 install .

.. tip::

   Unlike the console, the Python interface is fully compatible with ``mpi4py``.
   Both **OpenSn** and ``mpi4py`` share the same MPI communicator. Therefore,
   the Python module can be used in scripts that incorporate other tasks using
   ``mpi4py``.

Step 9 - Run Regression Tests
-----------------------------

To verify that the code compiled correctly, run the test scripts:

.. code-block:: shell

   cd /path/to/opensn
   test/run_tests -j<N>
   build/test/opensn-unit

Step 10 - OpenSn Documentation
------------------------------

If you configured the OpenSn build environment with support for building the
documentation (see **Step 8**), follow these instructions to install the
necessary tools and generate the documentation.

First, install the required Python packages using ``pip3``:

.. code-block:: shell

   pip3 install breathe myst-parser sphinx_rtd_theme

Then, from your ``build`` directory, generate the documentation with:

.. code-block:: shell

   cd build
   make doc

Once the build process is complete, you can view the generated documentation by
opening ``opensn/build/doc/index.html`` in your preferred web browser.
