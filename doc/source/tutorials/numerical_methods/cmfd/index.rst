CMFD Acceleration
=================

These tutorials show how to use coarse-mesh finite-difference (CMFD) acceleration for k-eigenvalue
problems and how its settings affect convergence and cost. Each one solves the C5G7 benchmark core
with unaccelerated power iteration or with several CMFD settings and compares the eigenvalue, the
number of power iterations and transport sweeps, and the wall time. Start with the basics tutorial.
The others can be read in any order.

The tutorials use a pin-homogenized model of the C5G7 core, in which every pin cell is one mesh
cell, so that each solve takes a few seconds. The benchmark tutorial solves the same core with
heterogeneous pins. Run each tutorial with four MPI processes and the ``pyopensn`` module. The
:doc:`CMFD user guide <../../../userguide/cmfd_acceleration>` describes every option in detail.
The method itself is fixed, so the tutorials vary only the coarse mesh, the coarse energy groups,
the transport update and the convergence tolerances.

The tutorials find the C5G7 cross sections by searching for ``test/assets/xs/c5g7`` in the
directory of the script, or the working directory of the notebook kernel, and in its parents. Keep
an exported script, or start the notebook kernel, inside an OpenSn checkout, for example in this
directory. Elsewhere the tutorial stops before it builds the mesh. To run a tutorial as a script
from this directory:

.. code-block:: shell

   jupyter nbconvert --to script cmfd_basics.ipynb
   mpiexec -n 4 python cmfd_basics.py

.. list-table::
   :header-rows: 1
   :widths: 30 70

   * - Tutorial
     - Topics
   * - :doc:`cmfd_basics`
     - Power iteration with and without CMFD, what one CMFD iteration does, status lines
   * - :doc:`cmfd_coarse_mesh`
     - ``coarse_mesh`` (``identity``, ``reactor_lattice``, ``cartesian``, ``logical_volumes``,
       ``local_aggregation``, ``global_aggregation``)
   * - :doc:`cmfd_energy_groupsets`
     - ``coarse_groups``, CMFD with more than one groupset
   * - :doc:`cmfd_transport_update`
     - ``sweeps_per_update``, Richardson and GMRES transport updates
   * - :doc:`cmfd_convergence_diagnostics`
     - ``k_tol`` and ``fission_source_tol``, ``verbose`` diagnostics, kept coarse values,
       skipped corrections, ``petsc_options``, adjoint problems
   * - :doc:`c5g7_benchmark`
     - The C5G7 core with heterogeneous pins and the pin fission-rate map

.. toctree::
   :maxdepth: 1

   cmfd_basics
   cmfd_coarse_mesh
   cmfd_energy_groupsets
   cmfd_transport_update
   cmfd_convergence_diagnostics
   c5g7_benchmark
