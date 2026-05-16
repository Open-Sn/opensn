=================
CMFD Acceleration
=================

Coarse-mesh finite difference (CMFD) acceleration can be used with
:py:class:`pyopensn.solver.PowerIterationKEigenSolver` through
:py:class:`pyopensn.solver.CMFDAcceleration`.

CMFD is a low-order scalar-flux acceleration method. The transport solve still
handles the angular flux, transport sweeps, material cross sections, and any
higher-order scattering moments. After each power-iteration transport update,
CMFD restricts the all-groups scalar flux and net face currents to a coarse
mesh, solves a coarse diffusion-like k-eigenvalue balance equation with
nonlinear current corrections, and prolongs a bounded scalar-flux correction
back to the transport solution.

Basic Usage
===========

Construct the transport problem as usual, then construct a CMFD acceleration
object and pass it to the power-iteration solver:

.. code-block:: python

   from pyopensn.solver import CMFDAcceleration, PowerIterationKEigenSolver

   cmfd = CMFDAcceleration(
       problem=phys,
       verbose=True,
   )

   solver = PowerIterationKEigenSolver(
       problem=phys,
       acceleration=cmfd,
       max_iters=400,
       k_tol=1.0e-7,
   )
   solver.Initialize()
   solver.Execute()

Supported Scope
===============

The CMFD accelerator currently supports:

* one or more groupsets,
* one shared all-groups CMFD coarse system applied after the AGS transport
  update,
* identity coarse meshes,
* rank-local same-material cell aggregation,
* reflecting and non-reflecting boundary treatment,
* scalar-flux correction limiting,
* P0 CMFD acceleration for transport problems with P0 or higher-order
  scattering.

The accelerator is traditional scalar CMFD. Higher-order transport scattering
is supported by the transport solve, but CMFD itself only accelerates the
zeroth moment scalar-flux balance. It does not accelerate higher angular
moments directly.

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

Important Parameters
====================

``coarse_mesh``
  Coarse-mesh construction method. Use ``"identity"`` for one CMFD coarse cell
  per transport cell or ``"local_aggregation"`` for rank-local aggregation.
  The default is ``"local_aggregation"``.

``aggregation_size``
  Target number of fine cells per local aggregated coarse cell. This only
  applies when ``coarse_mesh="local_aggregation"``. The default is ``16``.

``update_scheme``
  If ``True``, CMFD configures each groupset WGS solver as a loose transport
  update and applies a CMFD correction after each update. This is usually the
  fastest production mode for the current CMFD implementation.

``update_wgs_max_its`` and ``update_wgs_abs_tol``
  WGS controls used when ``update_scheme=True``. These are applied to every
  groupset.

``relaxation``
  Initial damping applied to the CMFD scalar-flux correction. The default is
  ``1.0``; adaptive relaxation may reduce or grow the starting value on later
  corrections.

``adaptive_relaxation``
  If ``True``, CMFD adjusts the next starting relaxation using recent accepted,
  damped, and rejected corrections. This is enabled by default so production
  inputs do not need to pick a single fixed damping value.

``adaptive_relaxation_min`` and ``adaptive_relaxation_max``
  Bounds for the adaptive starting relaxation.

``inactive_iterations``
  Number of initial power iterations before applying CMFD corrections. The
  transport update scheme is still active during these iterations. This warmup
  avoids applying a strong low-order correction to a poorly established
  high-order flux shape.

``correction_max_attempts``, ``correction_min_damping``, and
``negative_flux_tolerance``
  Correction-limiter controls. If a candidate correction produces non-finite
  fluxes, an invalid k-eigenvalue, or a scalar-flux undershoot beyond the
  allowed tolerance, CMFD halves the damping and retries. If no acceptable
  correction is found, the CMFD correction is skipped for that transport update.

``pi_max_its`` and ``pi_k_tol``
  Controls for the inner power iteration on the CMFD coarse system. The default
  ``pi_max_its=5`` is intended to give the coarse correction enough work to be
  useful without making the coarse solve dominate runtime.

``coarse_solver_policy``
  PETSc coarse-system solver policy. ``"auto"`` uses a direct PETSc LU solve
  below ``direct_coarse_solve_threshold`` global CMFD unknowns and GMRES with
  Jacobi above that threshold. ``"direct"``, ``"iterative"``, and
  ``"petsc_options"`` can be used to force a specific path.

``direct_coarse_solve_threshold``
  Maximum global CMFD unknown count for the direct PETSc LU path when
  ``coarse_solver_policy="auto"``. The default is ``5000``.

Diagnostics
===========

Set ``verbose=True`` on ``CMFDAcceleration`` to print stable ``CMFD_METRIC``
diagnostics. These are intended for regression tests and performance studies.

Coarse-mesh diagnostics include:

* ``global_fine_cells``
* ``global_coarse_cells``
* ``aggregation_ratio``
* ``global_unknowns``
* ``max_fine_cells_per_coarse_cell``
* ``max_faces_per_coarse_cell``
* ``undersized_coarse_cells``
* ``average_faces_per_coarse_cell``

Correction diagnostics include:

* ``damping``
* ``attempts``
* ``skipped``
* ``min_scalar_flux``
* ``nonfinite_flux``
* ``invalid_k``

Timing diagnostics include transport update time, restriction time, face-current
cache time, matrix assembly time, coarse solve time, correction update time, and
coarse KSP iteration counts.

Recommended Starting Point
==========================

For k-eigenvalue problems where normal power iteration is slow, start with:

.. code-block:: python

   cmfd = CMFDAcceleration(
       problem=phys,
       coarse_mesh="local_aggregation",
       aggregation_size=16,
       update_scheme=True,
       update_wgs_max_its=4,
       update_wgs_abs_tol=1.0e-4,
       relaxation=1.0,
       adaptive_relaxation=True,
       adaptive_relaxation_min=0.25,
       inactive_iterations=1,
       pi_max_its=5,
       coarse_solver_policy="auto",
   )

Then check:

* k-eigenvalue agreement against a converged non-accelerated or trusted
  reference solve,
* ``CMFD_METRIC c=correction m=skipped`` remains low or zero,
* ``nonfinite_flux`` and ``invalid_k`` remain zero,
* timing metrics show that reduced transport work pays for the coarse solve and
  coarse-mesh bookkeeping.

If convergence is noisy or corrections are skipped frequently, reduce
``relaxation`` or increase ``inactive_iterations``. If runtime is dominated by
transport sweeps, compare ``update_wgs_max_its`` values near the default. If
runtime is dominated by the coarse solve, compare ``coarse_solver_policy``
settings.
