Binding C++ classes to Python
=============================

The **OpenSn** package provides a set of Python bindings for its C++ classes.
The bindings are generated using the `pybind11 <https://github.com/pybind/pybind11>`_
library. Additional packages such as `NumPy <https://numpy.org>`_ are also
required to build the bindings.

Currently, there are two ways to access the Python bindings: either through a
**console** or through a **module**.

All binding code is located in the ``python/lib`` directory.

Namespaces
----------

In addition to the global namespace, **OpenSn** defines several internal
namespaces to organize related functionality:

- ``aquad`` — for angular quadrature
- ``diffusion`` — for diffusion solver
- ``fieldfunc`` — for field functions
- ``logvol`` — for logical volumes
- ``mesh`` — for mesh-related functionality
- ``math`` — for math-related functionality
- ``post`` — for post-processing
- ``response`` — for response evaluator
- ``settings`` — for changing program parameters (accessible only through the module)
- ``solver`` — for solvers
- ``source`` — for neutron point source and volumetric source
- ``xs`` — for cross sections

These namespaces should not be modified unless explicitly instructed by the
development team.

Classes
-------

Example of a binding of ``opensn::MyClass`` to the ``aquad`` namespace:

.. code-block:: cpp

   void
   WrapMyClass(py::module& aquad)
   {
     // clang-format off
     // my class
     auto my_class = py::class_<MyClass, std::shared_ptr<MyClass>, BaseClassIfBinded>(
       aquad,
       "MyClass",
       R"(
       Short description...

       Long description of MyClass...

       Wrapper of :cpp:class:`opensn::MyClass`.
       )"
     );
     my_class.def(
       py::init(
         [](int a, int b, const std::string & c)
         {
           return std::make_shared<MyClass>(a, b, c);
         }
       ),
       R"(
       Constructor documentation...

       Parameters
       ----------
       a: int
           ...
       b: int
           ...
       c: str
           ...
       )",
       py::arg("a"),
       py::arg("b"),
       py::arg("c")
     );
     my_class.def(
       "MyMethod1",
       &MyClass::MyMethod1,
       R"(
       Docstrings...

       Parameters
       ----------
       a: str
           String ...
       b: int, default=1
           Number of...
       )",
       py::arg("a"),
       py::arg("b") = 1
     );
     my_class.def_static(
       "MyClassMethod",
       &MyClass::MyClassMethod,
       R"(
       ...
       )",
       ...
     );
     // clang-format on
   }

   void
   py_aquad(py::module& pyopensn)
   {
     py::module aquad = pyopensn.def_submodule("aquad", "Angular quadrature module.");
     ...
     WrapMyClass(aquad);
   }

.. important::

   Docstrings must follow `Numpy docstring format <https://numpydoc.readthedocs.io/en/latest/format.html>`_.
   All non-empty lines in a docstring paragraph must have **the same indent level**.

.. tip::

   Math equations and figures can be added through the ``:math:`` and
   ``:img:`` directives.

Compilation for developers
--------------------------

While the Python module can be compiled and install using ``pip``, it is
recommended for developers to compile and install the module and the console
application inplace (i.e. in the source directory of **OpenSn**). This approach
will print all the CMake logs to the terminal, thus facilitating the debugging
process.

The compilation command for building the Python interface inplace is:

.. code-block:: shell

   python setup.py build_ext -i

Extra arguments can be appended to the build command, including:

* ``-h``, ``--help``: print the help message.
* ``-i``, ``--inplace``: ignore build-lib and put compiled extensions into the source directory alongside your pure Python modules.
* ``-g``, ``--debug``: compile in debug mode and turn off all optimizations.
* ``-f``, ``--force``: forcibly build everything (ignore file timestamps).
* ``-j<n>``, ``--parallel=<n>``: number of parallel build jobs. By default, it uses all threads of the builder processor.

Environment variables are used to set arguments for ``cmake`` builder:

* ``CMAKE_ARGS``: extra arguments to pass to CMake during configure step.
* ``BUILD_ARGS``: extra arguments to pass to CMake during build step.

.. code-block:: sh

   CMAKE_ARGS='-DOPENSN_WITH_DOCS=ON' BUILD_ARGS='--target __init__ --target doc' python setup.py build_ext --inplace

Pure Python plug-ins (future development)
-----------------------------------------

Pure Python plug-ins, such as submodule for plot utilities, can be added to the
``pyopensn`` folder. This feature enables seamless communication between the C++
API for intensive numerical computation and the Python API for flexibility and
data visualization.

.. important::

   This feature is only available through the module approach.
