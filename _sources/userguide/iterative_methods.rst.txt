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
- ``apply_tgdsa``
- ``tgdsa_l_abs_tol``
- ``tgdsa_l_max_its``
- ``tgdsa_verbose``
- ``tgdsa_petsc_options``

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
in the solution, together with a spectral-radius estimate.

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

``ags_convergence_check``
-------------------------

Allowed values are:

- ``"l2"``
- ``"pointwise"``

``l2`` is the default and is usually the simplest first choice.

``pointwise`` is stricter in a different sense: it tracks the maximum local
relative change in the flux vector. It can be useful when localized solution
changes matter more than the global norm.

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

In general, the groupset inner solves should be at least as tight as the AGS
solve, and usually somewhat tighter.

``max_ags_iterations``
----------------------

This is the hard cap on AGS passes.

What it tests:

- nothing by itself; it simply limits the total number of across-groupset
  iterations

Reasonable values:

- ``100`` is a good default
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
- ``tgdsa_l_abs_tol``
- ``tgdsa_l_max_its``
- ``tgdsa_verbose``
- ``tgdsa_petsc_options``

WGDSA
-----

WGDSA stands for within-group diffusion synthetic acceleration. Enabling
``apply_wgdsa`` creates a diffusion MIP solve over the current groupset and
applies the resulting correction to the updated flux.

Why use it:

- to accelerate scattering-dominated inner iterations
- especially useful with Richardson-style iteration
- can also help PETSc-based inner solves when used as part of the transport
  solve workflow

WGDSA controls:

- ``wgdsa_l_abs_tol``: stopping tolerance for the diffusion correction solve
- ``wgdsa_l_max_its``: hard cap on WGDSA iterations

Reasonable values:

- start with ``wgdsa_l_abs_tol=1.0e-4`` and ``wgdsa_l_max_its=30``
- if WGDSA is clearly helping but not converging tightly enough, reduce the
  tolerance to ``1.0e-5`` or raise the iteration cap modestly
- it is usually not necessary to solve the DSA system more tightly than the
  transport problem itself

TGDSA
-----

TGDSA stands for two-grid diffusion synthetic acceleration. Enabling
``apply_tgdsa`` constructs a collapsed one-group diffusion correction using
two-grid information derived from the groupset cross sections, then projects the
correction back to the groupset spectrum.

TGDSA controls:

- ``tgdsa_l_abs_tol``: stopping tolerance for the two-grid diffusion solve
- ``tgdsa_l_max_its``: hard cap on TGDSA iterations

Reasonable values:

- the same starting values as WGDSA are usually sensible:
  ``tgdsa_l_abs_tol=1.0e-4`` and ``tgdsa_l_max_its=30``
- TGDSA is a correction solve, so very tight tolerances are usually unnecessary

Why use it:

- to accelerate multigroup error modes that are not well treated by within-group
  acceleration alone

Practical guidance
------------------

For many problems:

- start without DSA
- if convergence is too slow, try WGDSA first
- add TGDSA when multigroup coupling remains difficult

DSA is most useful when there is an actual iteration problem to solve. It is
not automatically beneficial for every easy case.

Example:

.. code-block:: python

   {
       "groups_from_to": (0, 39),
       "angular_quadrature": quad,
       "inner_linear_method": "petsc_gmres",
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
- performs the transport solve through AGS and the groupset inner solves
- updates ``k_eff``
- repeats until the eigenvalue change is below ``k_tol`` or ``max_iters`` is
  reached

``k_tol`` is the stopping tolerance on the eigenvalue change between successive
power iterations. ``1.0e-8`` to ``1.0e-10`` is a reasonable production range,
with ``1.0e-10`` being common when a fairly tight ``k_eff`` is desired.

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
4. If multigroup convergence remains difficult, consider TGDSA and/or revisit
   groupset splitting.
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
