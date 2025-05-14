Installation guide
==================

These instructions provide setup guidance for both Linux (tested on Ubuntu
22.04.4 LTS) and macOS. While most steps are similar across platforms, some
adjustments may be needed depending on your specific OS version or
configuration.

.. note::

    - These instructions assume you are using the bash shell. If you use a
      different shell (e.g., zsh, fish), commands may require slight
      modifications.
    - Some steps require administrative privileges and the use of sudo. If you
      do not have such privileges, please consult your system administrator.

Install Development Tools
-------------------------

The following packages are required for **OpenSn** installation and development:

1. A recent version of ``clang++``/``g++`` that supports C++17
2. ``gfortran`` (required by the BLAS component of PETSc)
3. ``flex`` and ``bison`` (required by the PTSCOTCH component of PETSc)
4. Python3 v3.9+ and ``pip`` (required by PETSc and OpenSn)
5. Git version control system
6. CMake v3.12+
7. MPI (OpenMPI, MPICH, and MVAPICH have been tested)
8. HDF5
9. Doxygen (required for generating the OpenSn documentation)

.. tab-set::

   .. tab-item:: Linux (Ubuntu)

      .. code-block:: shell

         sudo apt install build-essential gfortran python3-dev python3-pip  \
             git cmake libopenmpi-dev flex bison libhdf5-mpi-dev doxygen
         export NPROC=$(nproc)

   .. tab-item:: macOS

      .. code-block:: shell

         brew install gcc python git cmake open-mpi flex bison hdf5-mpi doxygen
         export NPROC=$(sysctl -n hw.ncpu)

Install PETSc
-------------

The current supported version is `PETSc version 3.23.0+ <https://web.cels.anl.gov/projects/petsc/download/release-snapshots/petsc-3.23.0.tar.gz>`_.

In a scratch directory, install **PETSc** using the following commands:

.. code-block:: shell

   wget https://web.cels.anl.gov/projects/petsc/download/release-snapshots/petsc-3.23.0.tar.gz
   tar -xvzf petsc-3.23.0.tar.gz
   cd petsc-3.23.0

To configure **PETSc** with the required build options, run:

.. code-block:: shell

   ./configure  \
       --prefix=/path/to/dependencies/directory  \
       --download-hypre=1  \
       --with-ssl=0 --with-debugging=0  \
       --with-pic=1 --with-fc=0 --with-shared-libraries=1  \
       --download-f2cblaslapack=1  \
       --download-metis=1 --download-parmetis=1  \
       --download-superlu_dist=1  \
       --download-ptscotch=1  \
       --with-cxx-dialect=C++11  \
       --with-64-bit-indices \
       CC=mpicc CXX=mpicxx FC=mpifort \
       COPTFLAGS='-O3'  \
       CXXOPTFLAGS='-O3'  \
       FOPTFLAGS='-O3'  \
       PETSC_DIR=$PWD

If the configuration fails, consult **PETSc**'s user documentation.

Follow the **PETSc** build prompts to complete the **PETSc** installation.

After a successful install, export these environment variables every time a new
terminal is opened for launching applications with PETSc:

.. code-block:: shell

   export PETSC_ROOT="/path/to/dependencies/directory"
   export CMAKE_PREFIX_PATH="${PETSC_ROOT}/lib${CMAKE_PREFIX_PATH:+:${CMAKE_PREFIX_PATH}}"

Install VTK
-----------

Download `VTK 9.3.0 <https://www.vtk.org/files/release/9.3/VTK-9.3.0.tar.gz>`_ into a
suitable location.

Run the following to configure, build, and install:

.. code-block:: shell

   mkdir build && cd build
   cmake \
       -DCMAKE_INSTALL_PREFIX=/path/to/dependencies/directory  \
       -DBUILD_SHARED_LIBS=ON  \
       -DVTK_USE_MPI=ON  \
       -DVTK_GROUP_ENABLE_StandAlone=WANT  \
       -DVTK_GROUP_ENABLE_Rendering=DONT_WANT  \
       -DVTK_GROUP_ENABLE_Imaging=DONT_WANT  \
       -DVTK_GROUP_ENABLE_Web=DONT_WANT  \
       -DVTK_GROUP_ENABLE_Qt=DONT_WANT  \
       -DVTK_MODULE_USE_EXTERNAL_VTK_hdf5=ON  \
       -DCMAKE_BUILD_TYPE=Release  \
       ..
   make -j$NPROC
   make install

After a successful install, export these environment variables every time a new
terminal is opened for launching applications with VTK:

.. code-block:: shell

   export VTK_DIR="/path/to/dependencies/directory"
   export CMAKE_PREFIX_PATH="${VTK_DIR}/lib${CMAKE_PREFIX_PATH:+:${CMAKE_PREFIX_PATH}}"

Install Caliper
---------------

In a scratch directory, install **Caliper** using the following commands:

.. code-block:: shell

   wget https://github.com/LLNL/Caliper/archive/refs/tags/v2.10.0.tar.gz
   tar -xvzf v2.10.0.tar.gz
   cd Caliper-2.10.0
   mkdir build
   cd build
   cmake
       -DCMAKE_INSTALL_PREFIX=/path/to/dependencies/directory  \
       -DWITH_MPI=ON  \
       -DWITH_KOKKOS=OFF  \
       -DWITH_GOTCHA=OFF  \
       ..
   make -j$NPROC
   make install

After a successful install, export these environment variables every time a new
terminal is opened for launching applications with **Caliper**:

.. code-block:: shell

   export CALIPER_DIR="/path/to/dependencies/directory"
   export CMAKE_PREFIX_PATH="${CALIPER_DIR}/lib${CMAKE_PREFIX_PATH:+:${CMAKE_PREFIX_PATH}}"

Install Boost
-------------

In a scratch directory, install **Boost** using the following commands:

.. code-block:: shell

   wget https://archives.boost.io/release/1.86.0/source/boost_1_86_0.tar.gz
   tar -xvzf boost_1_86_0.tar.gz
   cd boost_1_86_0

Next, bootstrap the build system using:

.. code-block:: shell

   ./bootstrap.sh --prefix=/path/to/dependencies/directory --with-toolset=gcc

Build and install with:

.. code-block:: shell

   ./b2 install -j$NPROC --prefix=/path/to/dependencies/directory

After a successful install, export these environment variables every time a new
terminal is opened for launching applications with Boost:

.. code-block:: shell

   export Boost_ROOT=/path/to/dependencies/directory
   export LD_LIBRARY_PATH=$Boost_ROOT/lib:$LD_LIBRARY_PATH
   export C_INCLUDE_PATH=$BOOST_ROOT/include:$C_INCLUDE_PATH

Create virtual environment
--------------------------

It is recommended to install the Python interface of OpenSn within a virtual
environment rather than using the system-wide Python installation.

To create a virtual environment named ``opensn_env`` in
``/path/to/dependencies/directory``, execute:

.. code-block:: shell

   cd /path/to/dependencies/directory
   python3 -m venv opensn_env

To activate the virtual environment for the shell, execute:

.. code-block:: shell

   source /path/to/dependencies/directory/opensn_env/bin/activate

To deactivate the virtual environment, execute:

.. code-block:: shell

   deactivate

.. note::

   By default, Python searches for a package in the system-wide location first.
   In case of multiple installations of the same package, change ``sys.path``
   before importing.

Clone OpenSn
------------

.. important::

   If you want to contribute to **OpenSn**, it is strongly recommended that you
   first fork the **OpenSn** repository then clone your fork.

To clone the **OpenSn** repository:

.. code-block:: shell

   git clone https://github.com/Open-Sn/opensn.git

To clone your fork of **OpenSn**:

.. code-block:: shell

   git clone https://github.com/<username>/opensn.git

Configure and build OpenSn
--------------------------

**OpenSn** provides a Python interface. It is available in two formats: a
console application ``opensn`` and a Python module ``pyopensn``.

Classes and functions in the Python interface are detailed in :ref:`pyapi`.

.. attention::

   The console and the module are **not compatible** with each other. Attempting
   to import the module within the console will result in an import error. Users
   should select one approach and maintain consistent coding style throughout.

Console application
^^^^^^^^^^^^^^^^^^^

To compile the console application:

.. code-block:: shell

   mkdir build
   cd build
   cmake ..
   make -j$NPROC

.. danger::

   In the console application, all classes and functions are implicitly imported
   into the ``__main__`` module at startup. Therefore, omit submodule prefixes
   when referring to class or function names. Additionally, avoid redefining any
   **OpenSn** class or function names to prevent naming conflicts.

Module
^^^^^^

To compile the module and install in the Python ``site-packages`` path:

.. code-block:: bash

   pip install .

.. tip::

   Unlike the console, the Python interface is fully compatible with ``mpi4py``.
   Both **OpenSn** and ``mpi4py`` share the same MPI communicator. Therefore,
   the Python module can be used in scripts that incorporate other tasks using
   ``mpi4py``.

Run regression tests
--------------------

To verify that the implementation is fully compatible with the current API of
the code, run the test scripts:

.. code-block:: shell

   cd /path/to/opensn
   test/run_tests -d test/python -j$NPROC -v 1 -w 3

.. attention::

   Regression tests require both the console and the module. This can be
   achieved with:

   .. code-block:: shell

      cmake -DOPENSN_WITH_PYTHON_MODULE=ON ..

Build documentation
-------------------

Install the required Python packages to the virtual environment using ``pip``:

.. code-block:: shell

   cd doc
   pip install -r requirements.txt

.. important::

   Compiling documentation requires the Python module of OpenSn.

Then, from your ``build`` directory, generate the documentation with:

.. code-block:: shell

   cd build
   make doc

Once the build process is complete, you can view the generated documentation by
opening ``build/doc/html/index.html`` in your preferred web browser.
