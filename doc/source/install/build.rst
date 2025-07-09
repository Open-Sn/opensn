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

Install `pandoc <https://pandoc.org/installing.html>`_, and required Python
packages to the virtual environment using ``pip``:

.. code-block:: shell

   cd doc
   pip install -r doc/requirements.txt

.. important::

   Compiling documentation requires the **Python module of OpenSn**.

Then, from your ``build`` directory, generate the documentation with:

.. code-block:: shell

   make html

Once the build process is complete, you can view the generated documentation by
opening ``doc/build/html/index.html`` in your preferred web browser.
