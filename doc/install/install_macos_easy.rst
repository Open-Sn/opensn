Easy Install on macOS
=========================

The following instructions were tested on macOS Sonoma (version 14).

Step 1 - Install Development Tools
----------------------------------

The following packages should be installed via `Homebrew <https://brew.sh>`_:

- ``make``
- ``cmake``
- ``gcc``
- ``curl``
- ``bison``
- ``m4``
- ``doxygen``
- ``sphinx``

To install any missing packages, use the following command:

.. code-block:: shell

   brew install <package name>

Pay close attention to the output of ``brew``. It may be necessary to add
locations to your ``PATH`` environment variable in order to use the newly
installed utilities.

.. important::

   The ``sphinx`` package may not be available via ``Homebrew``. It is only
   required if you intend to build the OpenSn documentation.

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

MPICH recommends that the build and install directories be placed outside the
source tree (``/path/to/mpich``).

To check the installation, run:

.. code-block:: shell

   /path/to/mpich/install/bin/mpicc --version

If the installation was successful, a message similar to the following should
appear:

.. code-block:: text

   gcc-13 (Homebrew GCC 13.2.0) 13.2.0
   Copyright (C) 2023 Free Software Foundation, Inc.
   This is free software; see the source for copying conditions.  There is NO
   warranty; not even for MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.

Ensure that the message shows the desired compilers.

If the installation is successful, set the following environment variables in
your ``~/.bashrc`` or ``~/.bash_profile`` script:

.. code-block:: shell

   export MPI_DIR="/path/to/mpich/install"
   export PATH="${MPI_DIR}/bin:${PATH}"

   export CC="${MPI_DIR}/mpicc"
   export CXX="${MPI_DIR}/mpicxx"
   export FC="${MPI_DIR}/mpifort"
   export F77="${MPI_DIR}/mpif77"

Step 3 - Clone OpenSn
---------------------

.. note::

   If you want to contribute to **OpenSn**, it is strongly recommended to first
   fork the OpenSn repository into your own Git account and then clone your
   fork.

Clone the OpenSn repository or your fork:

.. code-block:: shell

   git clone https://github.com/Open-Sn/opensn.git /path/to/opensn

or

.. code-block:: shell

   git clone https://github.com/<username>/opensn.git /path/to/opensn

Step 4 - Set Up the Environment
-------------------------------

.. important::

   Xcode 15's linker introduces changes that may cause issues. It may be
   necessary to modify ``configure_dependencies.py`` to use PETSc 3.20.x. You
   may also need to add the following line to OpenSn's ``CMakeLists.txt``:

   ::

      link_libraries("-ld_classic")

Next, run the script to compile the necessary dependencies:

.. code-block:: shell

   mkdir -p /path/to/dependencies
   cd /path/to/opensn
   python3 tools/configure_dependencies.py -d /path/to/dependencies

It is recommended that ``/path/to/dependencies`` be located outside the OpenSn
source tree.

Set the environment variable for building OpenSn:

.. code-block:: shell

   export CMAKE_PREFIX_PATH=/path/to/dependencies${CMAKE_PREFIX_PATH:+:${CMAKE_PREFIX_PATH}}

Step 5 - Configure and Build OpenSn
-----------------------------------

Lua interface
^^^^^^^^^^^^^

OpenSn is configured within a build directory:

.. code-block:: shell

   cd /path/to/opensn
   mkdir build
   cd build
   cmake ..

To configure with support for building the documentation, use:

.. code-block:: shell

   cd /path/to/opensn
   mkdir build
   cd build
   cmake -DOPENSN_WITH_DOCS=ON ..

In general, the build directory will be located within the source tree.

Once configuration is complete, the Lua interface of OpenSn can be built within
the build directory using:

.. code-block:: shell

   make -j<N>

.. note::

   OpenSn may need to be reconfigured when dependencies change or new files are
   added. If this occurs, clear the ``build`` directory and repeat the
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

Step 6 - Run Regression Tests
-----------------------------

To check if the code compiled correctly, execute the test scripts:

.. code-block:: shell

   cd /path/to/opensn
   test/run_tests -j<N>
   build/test/opensn-unit

Step 7 - OpenSn Documentation
-----------------------------

If you configured the OpenSn build environment with support for building the
documentation (see **Step 5**), follow these steps to install the necessary
tools and generate the documentation.

First, install the required Python packages using ``pip3``:

.. code-block:: shell

   pip3 install breathe myst-parser sphinx_rtd_theme

Then, from your ``build`` directory, run the following command to generate the
documentation:

.. code-block:: shell

   cd build
   make doc

Once the build process has completed, you can view the generated documentation
by opening ``opensn/build/doc/index.html`` in your favorite web browser.
