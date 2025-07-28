OpenSn Installation
===================

Install Development Tools
-------------------------

The following packages are required for **OpenSn** installation and development:

1. A recent version of ``clang++``/``g++`` that supports C++20
2. ``flex`` (required by the PTSCOTCH component of PETSc)
3. Python3 v3.9+ and ``pip`` (required by PETSc and OpenSn)
4. Git version control system
5. CMake v3.20.2+
6. MPI (OpenMPI, MPICH, and MVAPICH have been tested)
7. PETSc 3.17.0+, Boost 1.86+, HDF5 1.14+, VTK 9.3.0+, Caliper 2.10+
8. pybind11 installed via ``pip``
9. Pandoc and Doxygen (required for generating the OpenSn documentation)

.. tab-set::

   .. tab-item:: Linux (Ubuntu)

      .. code-block:: shell

         sudo apt install build-essential python3-dev python3-pip  \
             git cmake libopenmpi-dev flex doxygen pandoc
         export NPROC=$(nproc)

   .. tab-item:: macOS

      .. code-block:: shell

         brew install gcc python git cmake open-mpi flex doxygen pandoc
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

For developers, it is recommended to use the following command to install the
additional packages required for running regression tests:

.. code-block:: bash

   pip install .[dev]

.. tip::

   Unlike the console, the Python interface is fully compatible with ``mpi4py``.
   Both **OpenSn** and ``mpi4py`` share the same MPI communicator. Therefore,
   the Python module can be used in scripts that incorporate other tasks using
   ``mpi4py``.

Run regression tests
--------------------

Run the regression tests to verify your installation:

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

   Compiling documentation requires the **Python module of OpenSn**.

Then, from your ``build`` directory, generate the documentation with:

.. code-block:: shell

   cd build
   make doc

Once the build process is complete, you can view the generated documentation by
opening ``build/doc/html/index.html`` in your preferred web browser.
