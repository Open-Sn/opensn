=================
CMFD Acceleration
=================

Coarse-mesh finite difference (CMFD) acceleration can be used with
:py:class:`pyopensn.solver.PowerIterationKEigenSolver` through
:py:class:`pyopensn.solver.CMFDAcceleration`.

CMFD is a low-order scalar-flux acceleration method. The transport solve still
handles angular fluxes, sweeps, material cross sections, and higher-order
scattering moments. After each power-iteration transport update, CMFD restricts
the scalar flux and face currents to a coarse system, solves a low-order
k-eigenvalue balance equation, and prolongs a bounded scalar-flux correction
back to the transport solution.

With CMFD, each power iteration performs a configured number of high-order WGS
transport update iterations, then applies a bounded low-order scalar-flux
correction. The defaults use one WGS update iteration, automatic current
closure, fixed correction relaxation, and a transport-current balance check
before outer power iteration is allowed to converge.

CMFD owns the WGS update controls through ``update_wgs_max_its`` and
``update_wgs_abs_tol``. The default is a one-iteration transport update before
each CMFD correction. To make each transport update more tightly converged
before CMFD is applied, increase ``update_wgs_max_its`` and tighten
``update_wgs_abs_tol``.

Basic Usage
===========

Construct the transport problem as usual, then pass a CMFD acceleration object
to the power-iteration solver:

.. code-block:: python

   from pyopensn.solver import CMFDAcceleration, PowerIterationKEigenSolver

   cmfd = CMFDAcceleration(problem=phys)

   solver = PowerIterationKEigenSolver(
       problem=phys,
       acceleration=cmfd,
       max_iters=400,
       k_tol=1.0e-6,
   )
   solver.Initialize()
   solver.Execute()

For large multigroup lattice problems, start by tuning only the spatial and
energy aggregation:

.. code-block:: python

   cmfd = CMFDAcceleration(
       problem=phys,
       current_closure="auto",
       aggregation_size=32,
       group_aggregation_size=(num_groups + 3) // 4,  # about 4 CMFD energy groups
       relaxation=0.5,
       update_wgs_max_its=1,
       update_wgs_abs_tol=1.0e-12,
       balance_residual_tolerance=1.0e-6,
   )

   solver = PowerIterationKEigenSolver(
       problem=phys,
       acceleration=cmfd,
       k_tol=1.0e-6,
       max_iters=400,
   )

The expression ``(num_groups + N - 1) // N`` chooses a
``group_aggregation_size`` that gives approximately ``N`` total CMFD energy
groups. For example, with 361 transport groups, ``(361 + 3) // 4`` gives a
group aggregation size of 91 and therefore 4 CMFD energy groups.

Common User Options
===================

These are the options most users should tune.

``coarse_mesh``
  Coarse-mesh construction method. The default is ``"local_aggregation"``.
  ``"local_aggregation"`` builds connected same-block aggregates locally on
  each MPI rank. ``"global_aggregation"`` builds connected same-block
  aggregates that may span MPI ranks and can reduce the coarse problem size on
  larger distributed meshes. ``"identity"`` creates one CMFD cell per transport
  cell and is mainly useful for debugging and method comparisons.

  Both aggregation modes are logical CMFD coarse-space constructions. They do
  not repartition the transport mesh, create a new mesh, or change
  the ownership of high-order transport cells. Transport sweeps, material
  lookup, and angular-flux communication continue to use the original fine-mesh
  partitioning.

  With ``"global_aggregation"``, every CMFD coarse cell has a single owning MPI
  rank even when its member fine cells span multiple ranks. Ownership is chosen
  from the rank containing the largest number of member fine cells, with ties
  resolved by the lower rank. Only the owner stores the full coarse cell and
  owns the corresponding PETSc rows. Ranks that own member fine cells store
  local fine-cell membership records so scalar fluxes, currents, and CMFD
  corrections can be restricted to and prolonged from the coarse-cell owner.

  ``"global_aggregation"`` currently has these restrictions:

  * a coarse cell contains only connected fine cells from the same mesh block;
  * connectivity is based on face-neighbor adjacency;
  * disconnected regions with the same block id become separate coarse cells;
  * ``aggregation_size`` is a target maximum, so boundary and disconnected
    aggregates may contain fewer fine cells;
  * no material homogenization across different blocks is performed.

  The global mode gathers enough distributed cell and face metadata to build
  the logical coarse graph consistently on all ranks. This is usually worthwhile
  when rank-local aggregation creates too many small coarse cells, but it can
  increase setup memory and communication compared with ``"local_aggregation"``.

``current_closure``
  Face-current closure used in the CMFD operator. Valid choices are ``"auto"``,
  ``"net"``, and ``"partial"``. The default is ``"auto"``. The automatic mode
  probes the early coarse-balance behavior and may choose net-current closure,
  partial-current closure, or a blend.

  ``"net"`` matches the signed transport current across each coarse face using
  a current-correction term in the coarse diffusion operator. ``"partial"``
  builds the face coupling from the outgoing partial currents on the two sides
  of the coarse face, normalized by the adjacent coarse scalar fluxes. The
  partial-current form can provide a stronger correction, but it is more
  sensitive to noisy or very small coarse fluxes. The automatic mode is the
  recommended starting point, but a fixed value is useful when comparing
  methods, reproducing a benchmark setting, or diagnosing a case where auto
  selection is not robust.

``aggregation_size``
  Target number of fine cells per aggregated CMFD coarse cell for
  ``"local_aggregation"`` and ``"global_aggregation"``. The default is ``32``.
  Larger values reduce the cost of the coarse solve but make the coarse
  correction less detailed. Smaller values increase coarse-system cost but can
  improve robustness. For lattice problems, values from roughly 8 to 64 are
  typical exploration points. Use wall time and final k agreement, not only the
  number of power iterations, when comparing values.

``group_aggregation_size``
  Number of transport energy groups per CMFD coarse energy group. This is not
  the final number of CMFD energy groups. The default is ``1``, which preserves
  the full transport group structure in the low-order system. To target a
  specific number of coarse groups, use::

     group_aggregation_size = (num_groups + target_coarse_groups - 1) // target_coarse_groups

  For example, ``group_aggregation_size=16`` means 16 transport groups per CMFD
  coarse group. By contrast, ``(num_groups + 15) // 16`` chooses an aggregation
  size that gives about 16 total CMFD coarse groups.

``relaxation``
  Requested relaxation factor for the scalar-flux correction. The default is
  ``0.5``. The correction limiter may reduce the damping for a single update or
  skip that update if the correction would produce non-finite values, an invalid
  k-eigenvalue, or an unacceptable scalar-flux undershoot. Smaller values are
  more conservative and can reduce oscillation; larger values can accelerate
  well-behaved cases but may trigger more damping or rejected corrections.

``update_wgs_max_its``
  Maximum WGS iterations used for each transport update. The default is ``1``.
  Increasing this can make each transport update more accurate, but also more
  expensive. The usual fast CMFD workflow keeps this small and lets CMFD provide
  the outer acceleration. For a tightly converged-transport workflow, increase
  this value and set a tight ``update_wgs_abs_tol``.

``update_wgs_abs_tol``
  WGS absolute tolerance used for each transport update. The default is
  ``1.0e-12``. When ``update_wgs_max_its=1``, this tolerance is normally not the
  stopping criterion because only one WGS iteration is allowed. It matters when
  ``update_wgs_max_its`` is large enough for the WGS solve to stop by residual
  tolerance before reaching the iteration limit.

``balance_residual_tolerance``
  Restricted transport-current balance residual required before CMFD permits the
  outer power iteration to converge. The default is ``1.0e-6``. This is an
  additional convergence check beyond ``k_eff_change`` and prevents false
  convergence when the eigenvalue appears stable but the CMFD-restricted
  transport balance is still inconsistent.

  This residual is a CMFD consistency guard, not an estimate of the
  k-eigenvalue error. It should not automatically be made tighter than the
  requested outer ``k_tol``. For large 3D problems, ``10*k_tol`` can be too
  strict and can force unnecessary iterations after the eigenvalue has already
  reached the requested accuracy. A useful starting point is to keep the
  default for production inputs, then loosen it only after comparing the final
  k-eigenvalue against a trusted PI or reference result.

Developer And Debug Options
===========================

The following options remain available for investigations, regression tests, and
method development. Production inputs should usually leave them at their
defaults.

``inactive_iterations``
  Number of initial power iterations before applying CMFD corrections. The
  default is ``0``. Transport update controls are still active during
  inactive iterations. Increase this only when debugging startup behavior where
  very early CMFD corrections are repeatedly rejected.

``correction_max_attempts``, ``correction_min_damping``, ``negative_flux_tolerance``
  Internal correction-limiter controls. These determine how many damping
  attempts are made and how much scalar-flux undershoot is tolerated before a
  correction is skipped. Leave these at their defaults unless you are debugging
  rejected corrections. Raising the allowed undershoot or lowering the minimum
  damping can hide instability; changing these should be accompanied by checks
  of final k and scalar-flux behavior.

``l_abs_tol`` and ``max_iters``
  Coarse linear solver tolerance and maximum iterations. Defaults are
  ``1.0e-7`` and ``100``. These control the CMFD coarse linear solves, not the
  outer transport convergence. CMFD corrections are not applied when a coarse
  linear solve fails to converge, because an unconverged low-order solve can
  produce a misleading correction. If this occurs with an iterative coarse
  solver, use a direct coarse solve for modest coarse systems, increase
  ``max_iters``, or provide stronger PETSc solver options. The default
  direct-solve threshold is intended to avoid automatic LU factorization of
  very large coarse systems; direct solves are limited by host memory and
  factorization cost.

``pi_max_its`` and ``pi_k_tol``
  Coarse k-eigenvalue power-iteration controls. Defaults are ``50`` and
  ``1.0e-8``. These control the inner CMFD coarse k solve and are separate from
  ``PowerIterationKEigenSolver(k_tol=...)``. If the coarse k solve is too loose,
  CMFD corrections can be noisy; if it is too tight, coarse-solve work can
  dominate runtime.

``coarse_solver_policy`` and ``direct_coarse_solve_threshold``
  PETSc coarse-solver selection controls. ``"auto"`` uses a direct PETSc LU
  solve for coarse systems below ``direct_coarse_solve_threshold`` global
  unknowns and an iterative GMRES/Jacobi solve otherwise. The default threshold
  is ``20000`` global unknowns. Force ``"direct"`` only for small or moderate
  coarse systems that fit comfortably in host memory. Use ``"iterative"`` when
  direct factorization is too expensive. If an iterative coarse KSP does not
  converge, CMFD skips the correction for that power iteration; increase
  ``max_iters``, use ``"petsc_options"`` with a stronger PETSc configuration,
  reduce coarse-system stiffness with finer aggregation, or switch to
  ``"direct"`` when memory permits.

``petsc_options``
  Additional PETSc options for developer experiments. This is normally empty.
  Options are intended for CMFD coarse-solver experiments, not for routine input
  tuning.

Multi-Groupset Behavior
=======================

CMFD is constructed as one accelerator for the transport problem, not one
accelerator per groupset. For problems with multiple groupsets, the transport
solver advances the groupsets through the normal AGS process, then CMFD
restricts the full all-groups scalar flux to one shared coarse system. The
coarse solve therefore couples all transport groups represented in the problem,
regardless of how those groups are partitioned into groupsets.

Diagnostics
===========

Normal outer-iteration verbosity prints compact PI and CMFD status lines. The
PI line reports the power-iteration state. When CMFD acceleration is active, a
separate line reports whether the CMFD-specific convergence gate is satisfied:

.. code-block:: text

   PI iteration = 24, k_eff = 0.9293662, k_eff_change = 5.74378e-09; WGS = iteration_limit
   CMFD convergence check = not_satisfied, transport-current balance residual = 2.300000e-06 (balance_residual_tolerance = 1.000000e-06)

The PI solver declares convergence only when both the eigenvalue tolerance and
the CMFD convergence check are satisfied. If a correction is rejected, the line
also reports ``correction = skipped`` with a reason, for example:

.. code-block:: text

   CMFD convergence check = not_satisfied, transport-current balance residual = 7.140531e-04 (balance_residual_tolerance = 1.000000e-05), correction = skipped (negative_flux_guard)

Set ``verbose=True`` on ``CMFDAcceleration`` to print stable ``CMFD_ACCEL``
diagnostic lines for early iterations and periodic later iterations. These are
intended for regression tests and performance studies.

Correction Safeties And Warnings
================================

CMFD never applies an untrusted coarse correction. If a safety check rejects a
correction, the current power iteration uses the unaccelerated transport update
instead. After repeated skipped corrections, CMFD returns the raw transport
k-eigenvalue update so that power iteration can continue to move forward even if
the low-order correction is temporarily unusable.

The most common skip reasons are:

``coarse_linear_solve_not_converged``
  One or more PETSc coarse linear solves did not converge within ``max_iters``.
  The correction is skipped because an unconverged low-order solve can move the
  high-order problem in the wrong direction. If the coarse system is modest,
  switch to ``coarse_solver_policy="direct"``. If direct is too expensive,
  increase ``max_iters`` or use ``coarse_solver_policy="petsc_options"`` with a
  stronger iterative PETSc configuration.

``negative_flux_guard``
  Every attempted correction, after damping, drove the scalar flux below the
  allowed undershoot. First reduce ``relaxation`` or use finer spatial/energy
  aggregation. Repeated negative-flux skips are usually a sign that the coarse
  correction is too aggressive for the current transport state. Changing
  ``negative_flux_tolerance`` is a developer/debug action and should be followed
  by checks of scalar flux and final k agreement.

``invalid_k`` or ``nonfinite_flux``
  The candidate correction produced a non-positive/non-finite k-eigenvalue or
  non-finite scalar flux. Treat this as an unstable correction: reduce
  ``relaxation``, use less aggressive aggregation, and check the coarse solver
  settings.

Warnings about skipped corrections are actionable. Occasional early skips can
be harmless because the unaccelerated transport update is retained. Persistent
skips usually mean that CMFD is not accelerating the calculation effectively and
the input should be adjusted.

Supported Scope
===============

The CMFD accelerator currently supports:

* one or more groupsets,
* one shared all-groups CMFD coarse system applied after the AGS transport
  update,
* identity coarse meshes,
* rank-local same-material cell aggregation,
* optional energy-group aggregation in the CMFD low-order system,
* reflecting and non-reflecting boundary treatment,
* scalar-flux correction limiting,
* P0 CMFD acceleration for transport problems with P0 or higher-order
  scattering.

The accelerator is scalar CMFD. Higher-order transport scattering is supported
by the transport solve, but CMFD itself accelerates only the zeroth moment
scalar-flux balance.

Limitations
===========

Current limitations are:

* CMFD is integrated with :py:class:`pyopensn.solver.PowerIterationKEigenSolver`,
  not :py:class:`pyopensn.solver.NonLinearKEigenSolver`.
* Local aggregation is rank-local. Coarse cells are not formed across MPI rank
  boundaries.
* The CMFD coarse solve uses PETSc on the host. GPU CMFD solves are not
  selected automatically.
* The correction is a scalar-flux correction. Angular flux and higher scattering
  moments remain governed by the high-order transport solve.
