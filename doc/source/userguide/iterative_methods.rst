=================
Iterative Methods
=================

Overview
========

OpenSn’s linear Boltzmann solver uses several levels of iteration. The most
important distinction for users is that not all iteration controls live in the
same place.

There are three levels to keep in mind:

1. inner groupset solves
2. across-groupset iterations
3. solver-specific outer iterations, such as k-eigenvalue iterations

In an OpenSn input, the inner groupset solver is selected in each groupset block,
while across-groupset controls are set in the problem ``options`` block.

Example:

.. code-block:: python

   phys = DiscreteOrdinatesProblem(
       mesh=grid,
       num_groups=20,
       groupsets=[
           {
               "groups_from_to": (0, 9),
               "angular_quadrature": quad,
               "inner_linear_method": "petsc_gmres",
               "l_abs_tol": 1.0e-6,
               "l_max_its": 200,
               "gmres_restart_interval": 30,
           },
           {
               "groups_from_to": (10, 19),
               "angular_quadrature": quad,
               "inner_linear_method": "petsc_gmres",
               "l_abs_tol": 1.0e-6,
               "l_max_its": 200,
               "gmres_restart_interval": 30,
           },
       ],
       xs_map=[{"block_ids": [0], "xs": xs}],
       options={
           "max_ags_iterations": 100,
           "ags_tolerance": 1.0e-6,
           "ags_convergence_check": "l2",
       },
   )

For most users, the practical questions are:

- which inner solver should be used in each groupset?
- should DSA be enabled?
- does the problem need multiple groupsets and AGS control?
- is the top-level solver a steady-state solve, a transient solve, a power
  iteration k-eigen solve, or a nonlinear k-eigen solve?


Where the Iteration Controls Live
=================================

Groupset controls
-----------------

These are specified inside each entry of ``groupsets=[...]``:

- ``inner_linear_method``
- ``l_abs_tol``
- ``l_max_its``
- ``gmres_restart_interval``
- ``allow_cycles``
- ``apply_wgdsa``
- ``wgdsa_l_abs_tol``
- ``wgdsa_l_max_its``
- ``wgdsa_verbose``
- ``wgdsa_petsc_options``
- ``wgdsa_solver_policy``
- ``wgdsa_direct_solve_threshold``
- ``wgdsa_boundary_closure``
- ``wgdsa_h_perp_convention``
- ``apply_tgdsa``
- ``tgdsa_l_abs_tol``
- ``tgdsa_l_max_its``
- ``tgdsa_verbose``
- ``tgdsa_petsc_options``
- ``tgdsa_solver_policy``
- ``tgdsa_direct_solve_threshold``
- ``tgdsa_boundary_closure``
- ``tgdsa_h_perp_convention``

These settings control the solve performed within a single groupset.

Problem-level controls
----------------------

These are specified in the problem ``options`` block:

- ``max_ags_iterations``
- ``ags_tolerance``
- ``ags_convergence_check``
- ``verbose_inner_iterations`` controls inner solver details, including WGS and AGS iterations.
- ``verbose_outer_iterations`` controls outer solver progress, including PI/NLKE iterations and
  transient step summaries.

These settings control how the solver coordinates multiple groupsets and how
much iteration information is printed.

A within-groupset solve that stops at its own iteration limit is treated as an
inexact inner solve: the across-groupset iteration continues until
``ags_tolerance`` is met or ``max_ags_iterations`` is reached. If
``max_ags_iterations`` is not set, each solver uses its own default: ``100``
for the steady-state and transient solvers, and ``1`` in each power iteration
of the power iteration k-eigen solver. A value set in the problem options is
used by every solver. The steady-state
and transient solvers print a warning when their transport solve ends without
converging (the across-groupset iteration with several groupsets, otherwise the
within-groupset solve).

Solver-level controls
---------------------

Some solvers add their own iteration settings.

For :py:class:`pyopensn.solver.PowerIterationKEigenSolver`:

- ``max_iters``
- ``k_tol``
- ``reset_phi0``

For :py:class:`pyopensn.solver.NonLinearKEigenSolver`:

- ``nl_abs_tol``
- ``nl_rel_tol``
- ``nl_sol_tol``
- ``nl_max_its``
- ``l_abs_tol``
- ``l_rel_tol``
- ``l_div_tol``
- ``l_max_its``
- ``l_gmres_restart_intvl``
- ``l_gmres_breakdown_tol``
- ``reset_phi0``
- ``num_initial_power_iterations``

Inner Groupset Methods
======================

The available groupset inner methods are:

- ``"classic_richardson"``
- ``"petsc_richardson"``
- ``"petsc_gmres"``
- ``"petsc_fgmres"``
- ``"petsc_bicgstab"``

These are selected with ``inner_linear_method`` in each groupset.

Example:

.. code-block:: python

   groupset = {
       "groups_from_to": (0, 9),
       "angular_quadrature": quad,
       "inner_linear_method": "petsc_gmres",
       "l_abs_tol": 1.0e-6,
       "l_max_its": 300,
       "gmres_restart_interval": 50,
   }

What these methods are solving
------------------------------

For a discrete ordinates groupset, the inner iteration repeatedly applies the
inverse transport sweep operator and updates the groupset flux moments until the
chosen convergence criterion is met.

The main difference between the available methods is how that fixed-point or
linear solve is managed:

- ``classic_richardson`` performs an explicit source iteration loop in OpenSn
  code
- the ``petsc_*`` methods use PETSc Krylov solvers wrapped around the sweep
  operator

As a practical rule, PETSc methods are usually the better starting point for
difficult problems. ``classic_richardson`` is simple and useful, but it is
generally the least robust option on strongly scattering or otherwise difficult
systems unless acceleration is also enabled.

``classic_richardson``
----------------------

This is OpenSn’s explicit source-iteration style method. On each iteration it:

- rebuilds the groupset source from the current iterate
- applies the inverse transport sweep
- optionally applies WGDSA and TGDSA corrections
- checks convergence using pointwise changes in ``phi`` and delayed angular flux

Why use it:

- simple algorithmic behavior
- useful for debugging and method studies
- pairs naturally with DSA for source-iteration style solves
- by far the least memory-intensive of the groupset inner solver methods

Why not use it by default:

- usually slower and less robust than Krylov methods on harder problems
- convergence can deteriorate badly for optically thick or highly scattering
  cases without acceleration
- in fixed-source and time-dependent solves, fission within the groupset is
  lagged by one iteration like scattering. WGDSA and TGDSA do not correct it, so
  a problem with a multiplication factor close to one converges at about that
  factor per iteration. The Krylov methods converge fission with the rest of the
  groupset operator.

Example:

.. code-block:: python

   {
       "groups_from_to": (0, 0),
       "angular_quadrature": quad,
       "inner_linear_method": "classic_richardson",
       "l_abs_tol": 1.0e-8,
       "l_max_its": 200,
   }

``petsc_richardson``
--------------------

This uses PETSc’s Richardson iteration around the transport operator.

Why use it:

- useful when a simple stationary iteration is desired but with PETSc-managed
  setup and convergence handling

A special case worth knowing:

- **if ``inner_linear_method="petsc_richardson"`` and ``l_max_its=1``, the code
  effectively performs a single sweep-style update**

For most problems, this is not the first PETSc method to try. GMRES is usually
the stronger default.

``petsc_gmres``
---------------

This uses PETSc GMRES around the transport operator.

For most steady-state transport problems, this is the best general-purpose
starting choice.

Why use it:

- usually the most robust of the groupset solver methods
- often the best choice for difficult scattering-dominated problems
- works with DSA preconditioning

Relevant control:

- ``gmres_restart_interval`` controls how many Krylov vectors GMRES keeps
  before it restarts

.. warning::

   GMRES can incur substantial memory usage because it stores restart vectors.
   Increasing ``gmres_restart_interval`` increases the Krylov storage footprint,
   so use larger values only when the convergence benefit justifies the extra
   memory.

Example:

.. code-block:: python

   {
       "groups_from_to": (0, 19),
       "angular_quadrature": quad,
       "inner_linear_method": "petsc_gmres",
       "l_abs_tol": 1.0e-6,
       "l_max_its": 300,
       "gmres_restart_interval": 30,
   }

``petsc_fgmres``
----------------

This uses PETSc flexible GMRES (FGMRES) around the transport operator, with
the preconditioner applied on the right.

Why use it:

- with ``apply_wgdsa`` or ``apply_tgdsa`` when the DSA diffusion system is
  solved iteratively. The iterative solve stops at the relative tolerance
  ``wgdsa_l_abs_tol``/``tgdsa_l_abs_tol``, so the DSA preconditioner differs
  slightly between applications. Standard GMRES assumes a fixed preconditioner
  and can then report convergence while the transport solution is still in
  error. Flexible GMRES remains exact. With standard GMRES, the transport error
  follows the DSA tolerance once that tolerance is looser than ``l_abs_tol``.
  A groupset that uses ``petsc_gmres`` or ``petsc_bicgstab`` with such a DSA
  tolerance prints a warning suggesting ``petsc_fgmres``. The warning compares
  only the tolerances. A direct diffusion solve, which the default ``"auto"``
  policy chooses for small systems on one MPI process, is the same in every
  application, and standard GMRES is then exact. Setting ``wgdsa_l_abs_tol`` and
  ``tgdsa_l_abs_tol`` no larger than ``l_abs_tol`` also avoids the error, at the
  cost of tighter diffusion solves.

It uses ``gmres_restart_interval`` like GMRES and stores about twice as many
restart vectors. Its residual test uses the unpreconditioned residual scaled by
the right-hand-side norm. Without a DSA preconditioner it behaves like
``petsc_gmres``.

``petsc_bicgstab``
------------------

This uses PETSc BiCGStab around the transport operator.

Why use it:

- as an alternative to GMRES when memory or restart behavior is undesirable
- as a solver-choice experiment for problems where GMRES is not ideal

In most common workflows, GMRES is still the more typical choice.


Convergence Controls
====================

Groupset tolerance and iteration count
--------------------------------------

Each groupset has:

- ``l_abs_tol``: default ``1.0e-6``
- ``l_max_its``: default ``200``

These control the inner solve stopping behavior.

For PETSc-based groupset methods, the residual printed in the log is the scaled
residual used by OpenSn’s custom convergence test for WGS solves.

For ``classic_richardson``, the convergence test is based on pointwise changes
in the solution, together with a spectral-radius estimate. The change of each
flux moment is measured relative to the scalar flux at the same node and group,
or to the moment itself when it is larger. The angular fluxes lagged at opposing
reflecting boundaries and broken sweep cycles must also change by less than
``l_abs_tol``, each measured relative to the scalar flux where the sweep reads
it, or to the lagged value itself when it is larger.

``l_abs_tol``
^^^^^^^^^^^^^

This is the main convergence tolerance for the groupset inner solve.

What it tests:

- for PETSc-based groupset methods, it tests the scaled WGS residual used by
  the OpenSn convergence check
- for ``classic_richardson``, it is the stopping tolerance for the pointwise
  flux-change test used by the Richardson iteration machinery

Reasonable values:

- ``1.0e-6`` is a good default for most production runs
- ``1.0e-4`` to ``1.0e-5`` is often acceptable for setup, debugging, or quick
  parameter studies
- ``1.0e-7`` to ``1.0e-8`` is reasonable when eigenvalues, localized
  responses, or balance quantities need tighter transport convergence

If the solve repeatedly reaches the outer iterations while the inner residual is
still large, the inner tolerance is probably too loose for that problem.

``l_max_its``
^^^^^^^^^^^^^

This is the hard cap on the number of groupset inner iterations.

What it tests:

- nothing by itself; it is a safeguard that stops the inner solve if the chosen
  convergence test is not met in time

Reasonable values:

- ``200`` is a good default
- ``100`` to ``300`` is a normal range for most problems
- higher values can be justified for very difficult scattering-dominated
  problems, but routinely hitting the cap usually means the method, groupset
  structure, or acceleration setup needs attention

What tolerance to choose
------------------------

As a practical rule:

- use looser tolerances while setting up and debugging a model
- tighten tolerances when k-eigenvalues, localized responses, or balance
  quantities need to be resolved accurately
- if outer iterations are still moving substantially, over-solving the inner
  iterations may not be the best use of work

This is especially important in multi-groupset and eigenvalue problems, where
the inner solve is only one part of the total iteration process.

GMRES restart interval
----------------------

``gmres_restart_interval`` only matters for ``petsc_gmres``.

Smaller restart values use less Krylov storage but may slow or degrade
convergence. Larger restart values can improve robustness but cost more memory
and work per restart cycle.

What it tests:

- it does not change the convergence test directly
- it controls how many GMRES Krylov vectors are retained before restart

Reasonable values:

- ``30`` is a good default
- ``20`` to ``50`` is a reasonable first range
- larger values such as ``60`` to ``100`` can help difficult problems if the
  extra memory is acceptable

.. note::

   On large production cases, the memory cost of GMRES restart vectors can be a
   practical limiter. If memory grows too much, reduce
   ``gmres_restart_interval`` first before assuming the only fix is more
   hardware.

``allow_cycles``
----------------

``allow_cycles`` controls whether sweep cycles are allowed in the groupset.

For most users, this is not the first iterative-method parameter to tune, but it
is part of the groupset solve setup and can matter on meshes whose sweep
dependencies contain cycles.

Across-Groupset Iteration
=========================

When a problem has multiple groupsets, OpenSn may perform across-groupset
iterations, usually abbreviated AGS.

The relevant problem options are:

- ``max_ags_iterations``
- ``ags_tolerance``
- ``ags_convergence_check``

What AGS does
-------------

AGS coordinates repeated solves of the groupsets until the full multigroup
solution is consistent across groupset boundaries.

In OpenSn, AGS:

- executes each WGS solver in sequence
- compares the new and previous full ``phi`` iterate
- uses either an ``l2`` or ``pointwise`` convergence check
- stops when the AGS tolerance is satisfied or the iteration limit is reached

Why this matters:

- if the problem has only one groupset, AGS is irrelevant
- if the problem has multiple groupsets, AGS is part of the main solve
  algorithm, not just a logging detail
- in fixed-source and time-dependent solves with fission, fission that emits
  from one groupset into another is iterated by AGS at about the multiplication
  factor of the problem per pass, which can take hundreds of passes (see
  :doc:`groupsets`)

``ags_convergence_check``
-------------------------

Allowed values are:

- ``"l2"``
- ``"pointwise"``

``l2`` is the default and is usually the simplest first choice. It measures
the L2 norm of the change of all flux moments over all spatial unknowns. The
norm is not normalized by the flux, so ``ags_tolerance`` is in the units of the
flux. Like the ``pointwise`` check, it stops when that change is below
:math:`(1-\varrho)\,` ``ags_tolerance``, where :math:`\varrho` is the estimated
convergence rate, so that the remaining error, about
:math:`\varrho/(1-\varrho)` times the last change, stays near the tolerance
even when the groupsets are strongly coupled.

``pointwise`` is stricter in a different sense: it tracks the maximum local
relative change in the flux vector, measuring each flux moment relative to the
local scalar flux or to the moment itself, whichever is larger. It can be useful
when localized solution changes matter more than the global norm.

.. _iterative_ags_tolerance:

``ags_tolerance``
-----------------

This is the stopping tolerance for the AGS iteration.

What it tests:

- the change between successive full multigroup ``phi`` iterates across AGS
  passes, using either the ``l2`` or ``pointwise`` check selected by
  ``ags_convergence_check``

Reasonable values:

- ``1.0e-6`` is a good default
- ``1.0e-5`` is often sufficient for debugging or coarse studies
- ``1.0e-7`` to ``1.0e-8`` may be appropriate when strong cross-groupset
  coupling must be converged tightly

In general, the groupset inner solves should be tighter than the AGS solve.
The change between AGS iterates includes the error that the inner solves
leave, so the AGS iteration cannot stop until that error is below its own test.
Compare the two in the same terms:

- With the PETSc inner methods, ``l_abs_tol`` bounds the relative residual of a
  groupset solve, not its error. A groupset whose iteration converges at rate
  :math:`\rho` leaves an error of up to about ``l_abs_tol`` :math:`/(1-\rho)`
  relative to the flux. :math:`\rho` is the convergence rate of the whole
  groupset iteration. Without acceleration it is set by everything the
  iteration lags: the within-group scattering, the scattering between the groups
  of the groupset and, in fixed-source and time-dependent solves, the fission
  within the groupset. Two groups that scatter strongly into each other converge
  slowly even without within-group scattering. WGDSA reduces the within-group
  part and TGDSA the part that scattering between the groups couples.
  ``classic_richardson`` already scales its test by :math:`1-\rho`.
- The ``l2`` check is absolute. ``ags_tolerance`` is in the units of the flux,
  so compare it with that relative error times the L2 norm of the flux. The
  ``pointwise`` check is relative.
- The AGS test also scales ``ags_tolerance`` by :math:`1-\varrho`.

One pass over the groupsets solves the coupling between them exactly when no
groupset receives a source from a groupset solved after it: no material scatters
from a later groupset into an earlier one (in an adjoint problem, with the
transposed scattering), no fission is emitted from one groupset into another,
and the problem does not use continuous slowing down. A power iteration holds
the fission source fixed during each transport solve, so there only the
scattering counts. The AGS iteration then stops after the first pass in which
every groupset solve converges, without testing the change between passes. The
accuracy of the result is then that of the groupset solves, so set ``l_abs_tol``
for the accuracy you need. A further pass would only continue the groupset
solves.

When the groupsets are coupled and the inner error is above the AGS test, the
AGS iteration does not stop at the inner tolerance. It runs slowly. A groupset
whose PETSc solve starts below ``l_abs_tol`` takes no Krylov iteration, but the
sweep that rebuilds its angular flux after the solve still applies one source
iteration to it. Each AGS pass
then advances that groupset by one source iteration, and the AGS iteration can
take tens of passes or reach ``max_ags_iterations``. The log shows this as
groupsets that converge at iteration 0 while the AGS change decreases slowly.
Tighten
``l_abs_tol``, add DSA to the groupsets with strong scattering, or loosen
``ags_tolerance``. Loosening ``ags_tolerance`` or switching to the
``pointwise`` check stops the iteration sooner, with the accuracy of the inner
solves.

``max_ags_iterations``
----------------------

This is the hard cap on AGS passes in one transport solve.

What it tests:

- nothing by itself; it simply limits the total number of across-groupset
  iterations

Defaults:

- if the option is not set, the steady-state and transient solvers use ``100``
  and the power iteration k-eigen solver uses ``1`` in each power iteration
- a value set in the problem options is used by every solver, including the
  power iteration solver, which then makes up to that many passes in every power
  iteration

Reasonable values for fixed-source and time-dependent solves:

- ``100``
- ``20`` to ``50`` is often enough for easier multi-groupset problems
- if this limit is reached regularly, the real issue is usually strong
  cross-groupset coupling, an overly aggressive groupset split, or the need for
  tighter inner solves

Example:

.. code-block:: python

   options = {
       "max_ags_iterations": 100,
       "ags_tolerance": 1.0e-6,
       "ags_convergence_check": "l2",
   }


DSA Options
===========

The groupset-level acceleration toggles are:

- ``apply_wgdsa``
- ``apply_tgdsa``

with associated controls:

- ``wgdsa_l_abs_tol``
- ``wgdsa_l_max_its``
- ``wgdsa_verbose``
- ``wgdsa_petsc_options``
- ``wgdsa_solver_policy``
- ``wgdsa_direct_solve_threshold``
- ``wgdsa_boundary_closure``
- ``wgdsa_h_perp_convention``
- ``tgdsa_l_abs_tol``
- ``tgdsa_l_max_its``
- ``tgdsa_verbose``
- ``tgdsa_petsc_options``
- ``tgdsa_solver_policy``
- ``tgdsa_direct_solve_threshold``
- ``tgdsa_boundary_closure``
- ``tgdsa_h_perp_convention``

WGDSA
-----

WGDSA stands for within-group diffusion synthetic acceleration. Enabling
``apply_wgdsa`` creates one MIP diffusion solve for the groupset. After each
groupset iteration, it corrects the scalar flux of every group with the
solution of a one-group diffusion problem whose source is the change of that
group's within-group scattering source. The groups are not coupled in the
correction. In a groupset with charged-particle groups and ``csda_enabled``,
WGDSA instead solves one diffusion problem per group in energy order, coupled
by in-groupset scattering and CSDA slowing-down (see :doc:`csda`). With
``classic_richardson`` the correction is applied after each sweep. With the
PETSc methods it is the preconditioner.

Why use it:

- to accelerate inner iterations that are slow because of within-group
  scattering in optically thick regions
- with any inner method: it is the correction step of Richardson iteration and
  the preconditioner of the Krylov methods

WGDSA controls:

- ``wgdsa_l_abs_tol``: residual tolerance of the iterative diffusion solve,
  relative to the norm of its right-hand side
- ``wgdsa_l_max_its``: maximum iterations of the iterative diffusion solve
- ``wgdsa_solver_policy``: PETSc solver setup policy for the WGDSA diffusion
  solve. ``"auto"`` uses direct PETSc LU on one MPI process for systems at or
  below ``wgdsa_direct_solve_threshold`` global unknowns, and conjugate
  gradients preconditioned by BoomerAMG otherwise. ``"direct"`` and
  ``"iterative"`` force those paths. ``"petsc_options"`` takes the KSP/PC setup
  from ``wgdsa_petsc_options``.
- ``wgdsa_direct_solve_threshold``: maximum global WGDSA diffusion unknown count
  for the automatic direct solve. The count is the number of spatial nodes times
  the number of groups in the groupset, or the number of spatial nodes for the
  per-group systems of a CSDA groupset. The default is ``20000``.
- ``wgdsa_petsc_options``: PETSc options for the diffusion solver. They must
  carry the solver's prefix (see :doc:`groupsets`).

The direct solve ignores ``wgdsa_l_abs_tol`` and ``wgdsa_l_max_its``.

Reasonable values:

- start with ``wgdsa_l_abs_tol=1.0e-4`` and ``wgdsa_l_max_its=30`` with
  ``classic_richardson`` or ``petsc_fgmres``
- with ``petsc_gmres`` or ``petsc_bicgstab`` and an iterative diffusion solve,
  use a ``wgdsa_l_abs_tol`` no larger than ``l_abs_tol`` (see
  ``petsc_fgmres`` above)
- if WGDSA is clearly helping but not converging tightly enough, reduce the
  tolerance to ``1.0e-5`` or raise the iteration cap modestly

TGDSA
-----

TGDSA stands for two-grid diffusion synthetic acceleration. Enabling
``apply_tgdsa`` builds, per material, the energy spectrum of the slowest error
mode of the groupset's iteration in energy and diffusion cross sections
collapsed with it. After each groupset iteration, it solves one one-group
diffusion problem whose source is the change of the scattering between the
groups of the groupset, and adds the solution times the spectrum to the scalar
flux of each group. It acts after WGDSA when both are enabled. TGDSA requires a
groupset with at least two groups and is not supported with CSDA.

TGDSA controls:

- ``tgdsa_l_abs_tol``: residual tolerance of the iterative two-grid diffusion
  solve, relative to the norm of its right-hand side
- ``tgdsa_l_max_its``: maximum iterations of the iterative two-grid diffusion
  solve
- ``tgdsa_solver_policy``: PETSc solver setup policy for the TGDSA diffusion
  solve, with the same choices and behavior as ``wgdsa_solver_policy``.
- ``tgdsa_direct_solve_threshold``: maximum global TGDSA diffusion unknown count
  for the automatic direct solve. The TGDSA system has one unknown per spatial
  node. The default is ``20000``.

Reasonable values:

- the same starting values as WGDSA are usually sensible:
  ``tgdsa_l_abs_tol=1.0e-4`` and ``tgdsa_l_max_its=30``, with the same caveat
  for ``petsc_gmres`` and ``petsc_bicgstab``
- TGDSA is a correction solve, so very tight tolerances are usually unnecessary

Why use it:

- to accelerate the error that scattering between the groups of a groupset
  couples across them, which within-group acceleration does not correct, for
  example in thermal groups with upscattering

Practical guidance
------------------

For many problems:

- start without DSA
- if convergence is too slow, enable WGDSA in groupsets whose groups scatter
  strongly within themselves
- enable TGDSA in a groupset that contains thermal upscattering, usually
  together with WGDSA
- keep every group that receives upscatter in the groupset that has TGDSA

See :ref:`groupsets_choosing_dsa` for when to use each correction and how to lay
out groupsets for TGDSA.

DSA is most useful when there is an actual iteration problem to solve. It is
not automatically beneficial for every easy case.

Example:

.. code-block:: python

   {
       "groups_from_to": (0, 39),
       "angular_quadrature": quad,
       "inner_linear_method": "petsc_fgmres",
       "l_abs_tol": 1.0e-6,
       "l_max_its": 200,
       "apply_wgdsa": True,
       "wgdsa_l_abs_tol": 1.0e-4,
       "wgdsa_l_max_its": 30,
       "wgdsa_verbose": False,
   }


Solver-Specific Outer Iteration
===============================

Steady-state source solves
--------------------------

For a standard steady-state source problem, the main iteration controls are the
groupset inner method and, when multiple groupsets are present, AGS.

Example:

.. code-block:: python

   ss_solver = SteadyStateSourceSolver(problem=phys)
   ss_solver.Initialize()
   ss_solver.Execute()

Transient solves
----------------

For transient problems, each timestep calls into the same underlying transport
iteration machinery. That means the groupset inner method and AGS settings still
matter, but now they matter per timestep.

If a transient run is unexpectedly expensive, the first iterative-method
settings to examine are usually:

- the chosen groupset inner method
- the inner tolerance
- whether angular flux storage is enabled for the transient formulation

Power iteration k-eigen solve
-----------------------------

:py:class:`pyopensn.solver.PowerIterationKEigenSolver` adds its own outer
iteration on top of the transport solve.

Relevant parameters:

- ``max_iters``
- ``k_tol``
- ``reset_phi0``

What it does:

- sets the fission source from the previous iterate
- performs the transport solve through AGS and the groupset inner solves, with
  one across-groupset pass unless the problem sets ``max_ags_iterations``
- updates ``k_eff``
- repeats until the relative eigenvalue change is below
  :math:`(1-\rho)\,` ``k_tol`` in two consecutive power iterations, or
  ``max_iters`` is reached

``k_tol`` is the tolerance on the estimated relative error of the eigenvalue.
A power iteration converging at rate :math:`\rho` changes :math:`k` by only
about :math:`1-\rho` times its error per iteration, so the relative change
:math:`|k_n - k_{n-1}| / k_n` must be below :math:`(1-\rho)\,` ``k_tol``, with
:math:`\rho` estimated from the change of the production-normalized flux. For a
problem with a dominance ratio close to one this takes more iterations than a
test on the change alone, which stops with an error of about :math:`1/(1-\rho)`
times ``k_tol``. The test must hold for two consecutive iterations, because an
accelerated eigenvalue sequence can overshoot and show one small change where it
turns around. ``1.0e-8`` to
``1.0e-10`` is a reasonable production range, with ``1.0e-10`` being common
when a fairly tight ``k_eff`` is desired.

With more than one groupset, each power iteration makes one across-groupset
pass unless the problem options set ``max_ags_iterations``, in which case it
makes up to that many. The fission source of a power iteration is itself only an
estimate, so converging the scattering coupling between groupsets inside each
power iteration does not change the converged eigenpair; it only adds sweeps.
With a single pass, the power iteration converges the coupling between
groupsets together with the fission source, and its rate estimate, taken from
the change of the flux, includes both. Strong upscattering split across many
groupsets can slow the power iteration; a few passes per power iteration then
reduce the number of power iterations.

``max_iters`` is the hard cap on power iterations. ``500`` to ``1000`` is a
reasonable starting range. If that cap is reached often, the issue is usually
not the cap itself but the overall transport convergence strategy.

Example:

.. code-block:: python

   keigen = PowerIterationKEigenSolver(
       problem=phys,
       max_iters=500,
       k_tol=1.0e-10,
   )
   keigen.Initialize()
   keigen.Execute()

Nonlinear k-eigen solve
-----------------------

:py:class:`pyopensn.solver.NonLinearKEigenSolver` uses a separate nonlinear
solve path.

Relevant parameters:

- nonlinear tolerances: ``nl_abs_tol``, ``nl_rel_tol``, ``nl_sol_tol``,
  ``nl_max_its``
- linear tolerances inside the nonlinear solve: ``l_abs_tol``, ``l_rel_tol``,
  ``l_div_tol``, ``l_max_its``
- GMRES controls for the nonlinear linearization solves:
  ``l_gmres_restart_intvl`` and ``l_gmres_breakdown_tol``
- optional warm start via ``num_initial_power_iterations``

Why use it:

- when the nonlinear eigenvalue solve is preferred over plain power iteration
- when convergence behavior of the nonlinear method is better suited to the
  problem

This is a different top-level algorithm from groupset WGS iteration. Its
internal linear solves are controlled by the ``l_*`` parameters on
``NonLinearKEigenSolver``, not by the groupset ``inner_linear_method`` or the
groupset ``l_*`` settings.

The nonlinear k-eigen solver supports multiple groupsets, but WGDSA/TGDSA
preconditioning is currently restricted to single-groupset nonlinear k-eigen
solves.

.. important::

   This is different from the other solver paths discussed in this section.
   ``NonLinearKEigenSolver`` uses its own internal linear solve. If a nonlinear
   k-eigen run needs tuning, start with the solver-level ``nl_*`` and ``l_*``
   parameters on :py:class:`pyopensn.solver.NonLinearKEigenSolver`, not with
   the groupset inner-solver settings.

What the main nonlinear controls do:

- ``nl_abs_tol``: absolute nonlinear residual tolerance
- ``nl_rel_tol``: relative nonlinear residual tolerance
- ``nl_sol_tol``: solution-update tolerance for the nonlinear solve
- ``nl_max_its``: hard cap on nonlinear iterations

Reasonable values:

- ``1.0e-8`` is a sensible starting point for ``nl_abs_tol`` and
  ``nl_rel_tol``
- ``25`` to ``50`` is a reasonable starting range for ``nl_max_its``
- ``nl_sol_tol`` is usually left at its default unless there is a specific
  reason to tune the solution-update criterion

What the internal linear controls do:

- ``l_abs_tol``: absolute tolerance for the linear solve inside each nonlinear
  step
- ``l_rel_tol``: relative tolerance for that same linear solve
- ``l_div_tol``: divergence guard for the linear solve
- ``l_max_its``: hard cap on linear iterations inside each nonlinear step
- ``l_gmres_restart_intvl``: restart size for the internal GMRES solve
- ``l_gmres_breakdown_tol``: breakdown guard for that GMRES solve

Reasonable values:

- ``1.0e-8`` is a reasonable starting value for ``l_abs_tol`` and
  ``l_rel_tol``
- ``50`` is a reasonable first value for ``l_max_its``
- ``30`` is a good default for ``l_gmres_restart_intvl``
- ``l_div_tol`` and ``l_gmres_breakdown_tol`` are usually left at their
  defaults unless the nonlinear linearization is clearly running into solver
  pathologies


Choosing a Method
=================

For most users, a sensible starting progression is:

1. Start with ``petsc_gmres``.
2. Use the default AGS settings unless there is a clear need to change them.
3. If convergence is slow on scattering-dominated problems, enable WGDSA.
4. If the groups of a groupset upscatter, enable TGDSA as well, and keep the
   whole upscatter range in that groupset (see :ref:`groupsets_choosing_dsa`).
5. Use ``classic_richardson`` when a simple source-iteration style method is
   specifically desired.

Keep in mind that ``petsc_gmres`` is not memory-free: the restart interval sets
how many Krylov vectors are retained between restarts, so aggressive values can
materially increase memory usage.

If the problem is easy, these choices may not matter much. If the problem is
difficult, the best improvements usually come from:

- a better inner method
- appropriate DSA
- sensible groupset decomposition
- tolerances that match the real accuracy target


Examples
========

Single groupset with GMRES
--------------------------

.. code-block:: python

   groupsets = [
       {
           "groups_from_to": (0, 9),
           "angular_quadrature": quad,
           "inner_linear_method": "petsc_gmres",
           "l_abs_tol": 1.0e-6,
           "l_max_its": 200,
           "gmres_restart_interval": 30,
       }
   ]

Classic Richardson with WGDSA
-----------------------------

.. code-block:: python

   groupsets = [
       {
           "groups_from_to": (0, 19),
           "angular_quadrature": quad,
           "inner_linear_method": "classic_richardson",
           "l_abs_tol": 1.0e-8,
           "l_max_its": 200,
           "apply_wgdsa": True,
           "wgdsa_l_abs_tol": 1.0e-4,
           "wgdsa_l_max_its": 30,
       }
   ]

Two groupsets with AGS control
------------------------------

.. code-block:: python

   phys = DiscreteOrdinatesProblem(
       mesh=grid,
       num_groups=40,
       groupsets=[
           {
               "groups_from_to": (0, 19),
               "angular_quadrature": quad,
               "inner_linear_method": "petsc_gmres",
               "l_abs_tol": 1.0e-6,
               "l_max_its": 200,
           },
           {
               "groups_from_to": (20, 39),
               "angular_quadrature": quad,
               "inner_linear_method": "petsc_gmres",
               "l_abs_tol": 1.0e-6,
               "l_max_its": 200,
           },
       ],
       xs_map=[{"block_ids": [0], "xs": xs}],
       options={
           "max_ags_iterations": 100,
           "ags_tolerance": 1.0e-6,
           "ags_convergence_check": "pointwise",
       },
   )

Power iteration k-eigen solve
-----------------------------

.. code-block:: python

   keigen = PowerIterationKEigenSolver(
       problem=phys,
       max_iters=1000,
       k_tol=1.0e-10,
       reset_phi0=True,
   )

Nonlinear k-eigen solve
-----------------------

.. code-block:: python

   keigen = NonLinearKEigenSolver(
       problem=phys,
       nl_abs_tol=1.0e-8,
       nl_rel_tol=1.0e-8,
       nl_max_its=50,
       l_abs_tol=1.0e-8,
       l_rel_tol=1.0e-8,
       l_max_its=50,
       l_gmres_restart_intvl=30,
       num_initial_power_iterations=5,
   )


Recommendations
===============

- Use ``petsc_gmres`` as the default first choice unless there is a specific
  reason to prefer another method.
- Remember that GMRES stores restart vectors, so aggressive restart sizes can
  significantly increase memory usage.
- Use ``classic_richardson`` when a simple source-iteration style solve is
  desired, especially for studies and debugging.
- Treat AGS settings as important only when the problem has multiple groupsets.
- Try WGDSA before TGDSA when the main issue is slow scattering convergence.
- Tighten tolerances only as far as the quantities of interest actually require.
- Revisit the iterative-method settings when changing problem difficulty,
  groupset structure, or solver type. A choice that is good for a fixed-source
  problem may not be the best one for a transient or k-eigen solve.
