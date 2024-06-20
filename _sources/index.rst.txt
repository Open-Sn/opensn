==================================
The OpenSn Discrete-Ordinates Code
==================================

OpenSn is an open-source, massively parallel, radiation transport code
designed to solve the discrete-ordinates Boltzmann transport equation
for neutral particles in a deterministic manner. OpenSn tackles a
variety of problems, including steady-state and dynamical
source-driven scenarios, as well as eigenvalue (criticality) problems.

Written in modern C++, OpenSn features a Lua interface, making it
versatile and accessible. From a simple laptop to a supercomputer,
OpenSn can be compiled and executed across a wide range of computing
platforms.

OpenSn aggregates decades of research and development in numerical
methods and parallel algorithms applied to radiation transport. The
most prominent features of OpenSn include:

-  A robust spatial finite-element method that is fully compatible with
   arbitrary polyhedral cells and adaptively refined grids.
-  Standard angular quadratures, as well as locally refined
   finite-element angular quadrature rules for the angular variable.
   Note that uncollided-flux treatment techniques to mitigate ray effects
   typically found in angular collocation methods will become available
   at a later date.
-  The traditional multigroup approximation applied to the energy
   variable. It is worth noting that multi-particle (e.g.,
   neutron+gamma) problems can be seamlessly handled.
-  Advanced iterative solution algorithms are employed to efficiently
   solve the resulting system of equations, including transport-sweep
   preconditioning with various synthetic accelerators.
-  Execution of massively parallel simulations on arbitrarily
   partitioned computational domains, utilizing efficient data
   aggregation and workflows to maximize performance. The code scales
   efficiently on a large number of MPI ranks, supporting up to 100,000+
   processes.

.. admonition:: Recommended publication for citing
   :class: tip

   Auhtors, "`OpenSn: title <https://doi.org/>`_,"
   *journal*, **volume**, page--page (year).

.. only:: html

   --------
   Contents
   --------

.. toctree::
   :maxdepth: 2

   getting_started
   theory/index
   tutorials/index
   coding_standard.md
   workflow.md
