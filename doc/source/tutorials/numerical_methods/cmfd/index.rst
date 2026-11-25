CMFD Acceleration
=================

These tutorials show how to use coarse-mesh finite-difference (CMFD) acceleration for k-eigenvalue
problems and how its options affect convergence and cost. Each one solves the C5G7 benchmark core
with unaccelerated power iteration or with several CMFD settings and compares the eigenvalue, the
number of power iterations and transport sweeps, and the wall time. Start with the basics tutorial.
The others can be read in any order.

The tutorials use a pin-homogenized model of the C5G7 core, in which every pin cell is one mesh
cell, so that each solve takes a few seconds. The benchmark tutorial solves the same core with
heterogeneous pins. Run each tutorial with four MPI processes and the ``pyopensn`` module. The
:doc:`CMFD user guide <../../../userguide/cmfd_acceleration>` describes every option in detail.

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
     - ``coarse_mesh`` (``identity``, ``local_aggregation``, ``global_aggregation``),
       ``aggregation_size``, ``allow_mixed_material_aggregation``
   * - :doc:`cmfd_energy_groupsets`
     - ``group_aggregation_size``, CMFD with more than one groupset
   * - :doc:`cmfd_closure_relaxation`
     - ``current_closure`` (``partial`` and ``net``) on thin and thick coarse cells,
       ``relaxation``, correction safeties
   * - :doc:`cmfd_transport_update`
     - ``update_wgs_max_its``, ``update_wgs_abs_tol``, Richardson and GMRES transport updates
   * - :doc:`cmfd_convergence_diagnostics`
     - ``balance_residual_tolerance``, the coarse eigenvalue solve (``pi_max_its``,
       ``pi_k_tol``, ``pi_phi_tol``), ``coarse_solver_policy``, ``verbose`` diagnostics, skipped
       corrections, ``inactive_iterations``, adjoint problems
   * - :doc:`c5g7_benchmark`
     - The C5G7 core with heterogeneous pins and the pin fission-rate map

.. toctree::
   :maxdepth: 1

   cmfd_basics
   cmfd_coarse_mesh
   cmfd_energy_groupsets
   cmfd_closure_relaxation
   cmfd_transport_update
   cmfd_convergence_diagnostics
   c5g7_benchmark
