Easy Install using built-in CMake script
========================================

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

Install dependencies
--------------------

Dependencies can be installed to ``/path/to/dependencies/directory`` with CMake
as follows:

.. code-block:: shell

   mkdir build_deps && cd build_deps
   cmake -DCMAKE_INSTALL_PREFIX=/path/to/dependencies/directory  \
       /path/to/opensn/tools/dependencies
   cd ..
   rm -rf build_deps

.. include:: build.rst
