Solvers
=======

The OpenSn solver interface is built around two components:

* a *problem* object that defines the transport model, mesh, groupsets,
  materials, sources, and boundary conditions, and
* a *solver* object that decides how that problem is advanced or iterated.

For linear Boltzmann transport, the most important problem classes are:

* :py:class:`pyopensn.solver.DiscreteOrdinatesProblem`
* :py:class:`pyopensn.solver.DiscreteOrdinatesCurvilinearProblem`

The main solver classes are:

* :py:class:`pyopensn.solver.SteadyStateSourceSolver`
* :py:class:`pyopensn.solver.TransientSolver`
* :py:class:`pyopensn.solver.PowerIterationKEigenSolver`
* :py:class:`pyopensn.solver.NonLinearKEigenSolver`

This section describes what each solver does, how the problem and solver layers
work together, which options belong where, and how transient problems can be
run either with :py:meth:`Execute` or with an explicit Python timestep loop.

.. note::

   GPU support is chosen on the problem side with ``use_gpus=True`` together
   with a compatible sweep path. In practice, solver GPU support depends on
   whether the associated problem configuration supports GPU sweeps.

Overview
========

The usual workflow is:

1. Construct a problem.
2. Construct a solver using that problem.
3. Call :py:meth:`Initialize`.
4. Call :py:meth:`Execute`, or for transient problems, repeatedly call
   :py:meth:`Advance`.

Example:

.. code-block:: python

   from pyopensn.solver import DiscreteOrdinatesProblem, SteadyStateSourceSolver

   phys = DiscreteOrdinatesProblem(...)

   solver = SteadyStateSourceSolver(problem=phys)
   solver.Initialize()
   solver.Execute()

.. note::

   The problem owns the transport state. The solver owns the algorithm used to
   advance that state. This distinction matters because you can often reuse the
   same problem with different solvers or switch the problem between
   steady-state, time-dependent, or adjoint modes.

.. note::

   In practice, most user scripts spend more lines constructing the problem than
   constructing the solver. That is expected. The problem captures the physics
   model; the solver choice mainly determines how that model is driven.

Base Classes
============

All solver objects derive from :py:class:`pyopensn.solver.Solver`, which
provides:

* :py:meth:`Initialize`
* :py:meth:`Execute`
* :py:meth:`Advance`

All transport problems derive from :py:class:`pyopensn.solver.Problem`, and all
linear Boltzmann problems derive from :py:class:`pyopensn.solver.LBSProblem`.

In practice, users usually work with the concrete classes rather than the base
classes, but it is useful to remember that:

* problem methods operate on the transport model and solution state
* solver methods operate on how that state is computed

.. note::

   If something sounds like "what problem am I solving?", it probably belongs on
   the problem object. If it sounds like "how do I advance or iterate that
   problem?", it probably belongs on the solver object.

.. note::

   This split is especially helpful when debugging. If results look physically
   wrong, inspect the problem configuration first. If results look physically
   reasonable but convergence or runtime is poor, inspect the solver and
   groupset iteration settings next.

LBS Problem
===========

:py:class:`pyopensn.solver.LBSProblem` is the shared base for the linear
Boltzmann problem classes. The Python API on this base class provides common
operations such as:

* scalar flux field-function access with
  :py:meth:`pyopensn.solver.LBSProblem.GetScalarFluxFieldFunction`
* derived field-function access with
  :py:meth:`pyopensn.solver.LBSProblem.CreateFieldFunction`
* current time and timestep queries with :py:meth:`GetTime` and
  :py:meth:`GetTimeStep`
* fission diagnostics with :py:meth:`ComputeFissionRate` and
  :py:meth:`ComputeFissionProduction`
* direct access to local scalar-flux vectors with :py:meth:`GetPhiOldLocal` and
  :py:meth:`GetPhiNewLocal`
* flux/source I/O helpers
* source and material reassignment with :py:meth:`SetPointSources`,
  :py:meth:`SetVolumetricSources`, and :py:meth:`SetXSMap`
* transport-mode changes with :py:meth:`SetAdjoint` and :py:meth:`SetForward`

.. note::

   These common methods are available regardless of which specific LBS solver
   you pair with the problem. This is why postprocessing and source/material
   updates are usually described at the problem level rather than at the solver
   level.

Discrete Ordinates Problem
==========================

:py:class:`pyopensn.solver.DiscreteOrdinatesProblem` is the standard Cartesian
discrete ordinates problem class.

Its constructor takes the main transport-model inputs:

* ``mesh``
* ``num_groups``
* ``groupsets``
* ``xs_map``
* ``boundary_conditions``
* ``point_sources``
* ``volumetric_sources``
* ``options``
* ``sweep_type``
* ``time_dependent``
* ``use_gpus``

The most important constructor inputs are:

* ``groupsets``: defines the angular quadrature and iterative behavior for
  ranges of energy groups
* ``xs_map``: maps mesh block ids to cross sections
* ``boundary_conditions``: defines vacuum, reflecting, isotropic, or arbitrary
  inflow boundaries
* ``options``: controls restart behavior, verbosity, precursor usage, angular
  flux storage, AGS behavior, and other global settings

Problem options available in Python include:

* ``max_mpi_message_size``
* ``restart_writes_enabled``
* ``write_delayed_psi_to_restart``
* ``read_restart_path``
* ``write_restart_path``
* ``write_restart_time_interval``
* ``use_precursors``
* ``use_source_moments``
* ``save_angular_flux``
* ``adjoint``
* ``verbose_inner_iterations``
* ``verbose_outer_iterations``
* ``max_ags_iterations``
* ``ags_tolerance``
* ``ags_convergence_check``
* ``verbose_ags_iterations``
* ``field_function_prefix_option``
* ``field_function_prefix``

``use_precursors`` controls whether delayed-neutron precursor treatment is kept
active for the problem. The default is ``True``. This should usually stay
enabled for transient and k-eigen workflows unless you explicitly want a
prompt-only model.

The setting is treated as user intent and persists across later
:py:meth:`SetXSMap` calls, even if the current cross-section map temporarily has
no precursor-bearing material. If cross sections are swapped, existing
precursor concentrations are remapped by local cell and precursor-family index;
new families start at zero and removed families are discarded. If a cell passes
through a material with zero precursors, its precursor history is dropped and
any later reintroduced precursor families restart from zero.

If any fissionable material in the active map contains precursor data and
``use_precursors=True``, then all fissionable materials in that map must
contain precursor data. Non-fissionable materials may have zero precursors.

.. note::

   Restart state is primarily problem-owned. The common restart payload stores
   flux moments, time metadata, precursor data, and any required
   discrete-ordinates angular state. Solver types can add only the small extra
   pieces they need on top of that common payload. Restart reads occur during
   :py:meth:`Initialize`. Timed restart dumps can be written by solvers that
   execute over multiple outer iterations or timesteps, such as
   :py:class:`pyopensn.solver.PowerIterationKEigenSolver` and
   :py:class:`pyopensn.solver.TransientSolver`.

.. note::

   Most numerical iteration settings that control transport sweeps live in the
   ``groupsets`` entries, not on the outer solver object. The main exception is
   :py:class:`pyopensn.solver.NonLinearKEigenSolver`, which has its own internal
   nonlinear and linear tolerance settings.

.. note::

   A good first mental split is:

   * ``groupsets`` control transport iteration and angular treatment,
   * ``xs_map`` and sources define the physics being solved,
   * ``options`` control global problem behavior such as restarts, verbosity,
     precursor treatment, and stored output.

Time-Dependent and Steady-State Modes
-------------------------------------

:py:class:`pyopensn.solver.DiscreteOrdinatesProblem` can operate in either
steady-state or time-dependent mode.

Use:

* :py:meth:`SetTimeDependentMode`
* :py:meth:`SetSteadyStateMode`
* :py:meth:`IsTimeDependent`

Important behavior:

* transient solves require the problem to be in time-dependent mode
* steady-state and k-eigenvalue solves require the problem to be in steady-state
  mode
* switching modes preserves sources and boundary conditions, unlike
  :py:meth:`SetAdjoint`, which performs a destructive mode reset

.. note::

   If you are moving from a steady-state solve into a transient solve, a common
   pattern is:

   1. solve the steady-state problem,
   2. call :py:meth:`SetTimeDependentMode`,
   3. construct and initialize a :py:class:`pyopensn.solver.TransientSolver`.

.. note::

   Mode is part of the problem state, not part of the solver object. That is
   why a solver may reject a problem that is otherwise fully defined if the
   problem is currently in the wrong mode.

.. note::

   For time-dependent mode, the problem must have been created with
   ``save_angular_flux=True``.

Angular-Flux-Specific Features
------------------------------

:py:class:`pyopensn.solver.DiscreteOrdinatesProblem` also provides:

* :py:meth:`SetBoundaryOptions`
* :py:meth:`GetPsi`
* :py:meth:`GetAngularFieldFunctionList`
* :py:meth:`ComputeLeakage`

Important requirement:

* if you need stored angular fluxes or angular flux field functions, enable
  ``options={"save_angular_flux": True}`` when creating the problem

.. note::

   ``save_angular_flux=True`` is primarily a capability switch. It enables
   workflows such as transient solves, angular-field output, and leakage
   postprocessing, but it also increases memory use.

DiscreteOrdinatesCurvilinearProblem
-----------------------------------

:py:class:`pyopensn.solver.DiscreteOrdinatesCurvilinearProblem` is the
curvilinear companion to the Cartesian problem class.

It uses the same general structure as
:py:class:`pyopensn.solver.DiscreteOrdinatesProblem`, but is intended for
curvilinear geometry and currently requires:

* ``coord_system=2`` for cylindrical coordinates

This problem type is experimental, and it should be treated that way in
production workflows.

.. note::

   The curvilinear problem class is best used when the mesh and quadrature are
   already being chosen with cylindrical geometry in mind. It is not just a
   cosmetic variant of the Cartesian problem.

Steady-State Source Solver
==========================

:py:class:`pyopensn.solver.SteadyStateSourceSolver` is the standard source-driven
steady-state transport solver.

Use it when:

* the problem is not time-dependent, and
* you want the solution driven by fixed volumetric sources, point sources, and
  boundary conditions rather than an eigenvalue solve

Constructor
-----------

The constructor takes:

* ``problem``: an existing :py:class:`pyopensn.solver.DiscreteOrdinatesProblem`

GPU support
-----------

This solver can run with GPU acceleration when the associated
:py:class:`pyopensn.solver.DiscreteOrdinatesProblem` is configured for GPU
sweeps. In practice, that means:

* ``use_gpus=True`` on the problem,
* ``sweep_type="AAH"``, and
* a non-curvilinear, non-time-dependent problem.

Example:

.. code-block:: python

   phys = DiscreteOrdinatesProblem(...)
   solver = SteadyStateSourceSolver(problem=phys)
   solver.Initialize()
   solver.Execute()

How It Operates
---------------

At a high level, the steady-state solver:

1. checks that the problem is in steady-state mode,
2. optionally reads restart data,
3. calls the AGS solver attached to the problem,
4. optionally writes restart data,
5. optionally computes precursor fields,
6. reorients the solution if adjoint mode is active, and
7. updates the completed transport state and balance information

The corresponding balance summary is available from
:py:meth:`ComputeBalanceTable`.

.. note::

   This is the simplest outer solver in the LBS stack. Most of the numerical
   work is still being done by the groupset and AGS machinery configured on the
   problem.

.. note::

   If a steady-state solve is unexpectedly slow, the most likely tuning targets
   are still the groupset iterative settings and acceleration options on the
   problem, not the steady-state solver object itself.

.. note::

   This solver is the standard baseline for source-driven transport. If the problem
   is not intrinsically an eigenvalue or transient problem, this is usually the
   first solver to reach for.

Adjoint Steady-State Solves
---------------------------

The same steady-state solver is also used for adjoint transport solves.

An adjoint solution answers a different question than a forward solution.
Instead of computing the particle field produced by a physical source, it
computes the importance of particles to a response of interest. That is why
adjoint solves are useful in detector analysis, response calculations,
importance studies, and related sensitivity workflows.

In OpenSn, an adjoint solve is still a steady-state transport solve, but the
problem is placed in adjoint mode so that the transport operator, materials,
and source interpretation are all treated in the adjoint sense.

There are two common ways to enable adjoint mode:

* set ``options={"adjoint": True}`` when constructing the problem, or
* call :py:meth:`pyopensn.solver.LBSProblem.SetAdjoint(True)` before solving

Typical pattern:

.. code-block:: python

   phys = DiscreteOrdinatesProblem(
       ...,
       options={"adjoint": True},
   )

   solver = SteadyStateSourceSolver(problem=phys)
   solver.Initialize()
   solver.Execute()

or:

.. code-block:: python

   phys = DiscreteOrdinatesProblem(...)
   phys.SetAdjoint(True)
   phys.SetVolumetricSources(volumetric_sources=adjoint_sources)
   phys.SetBoundaryOptions(boundary_conditions=adjoint_boundaries)

   solver = SteadyStateSourceSolver(problem=phys)
   solver.Initialize()
   solver.Execute()

Important behavior:

* adjoint mode is a property of the problem, not of the solver object
* if :py:meth:`SetAdjoint(True)` changes the current mode, OpenSn performs a
  destructive mode reset
* that reset reinitializes materials in adjoint mode, clears sources and
  boundary conditions, and zeroes scalar and angular flux state
* after such a mode change, the desired adjoint sources and boundary conditions
  must be reapplied before solving

At the end of the steady-state solve, OpenSn reorients the adjoint solution for
the discrete-ordinates problem so the stored result is in the expected adjoint
form.

.. note::

   Time-dependent adjoint problems are not supported. Adjoint mode is therefore
   a steady-state capability in the current solver stack.

.. note::

   A useful mental model is: the solver algorithm does not become a different
   object in adjoint mode. Instead, the same steady-state solve is applied to a
   problem whose transport data and source interpretation have been switched to
   the adjoint form.

Transient Solver
================

:py:class:`pyopensn.solver.TransientSolver` advances a
:py:class:`pyopensn.solver.DiscreteOrdinatesProblem` in time.

Use it when:

* the problem is in time-dependent mode, and
* you want backward Euler, Crank-Nicolson, or another theta-scheme transient
  update

GPU support
-----------

Transient solves do not currently support GPU acceleration.

Constructor
-----------

The transient solver constructor takes:

* ``problem``: an existing
  :py:class:`pyopensn.solver.DiscreteOrdinatesProblem`
* ``dt``: timestep size, default ``2.0e-3``
* ``stop_time``: final execution time, default ``0.1``
* ``theta``: time differencing parameter, default ``0.5``
* ``initial_state``: ``"existing"`` or ``"zero"``
* ``verbose``: enable or disable transient log output

.. important::

   Transient problems require ``save_angular_flux=True`` in the problem
   options. Set this on the
   :py:class:`pyopensn.solver.DiscreteOrdinatesProblem` before constructing the
   :py:class:`pyopensn.solver.TransientSolver`.

The module also exposes:

* ``pyopensn.solver.BackwardEuler = 1.0``
* ``pyopensn.solver.CrankNicolson = 0.5``

Meaning of ``initial_state``:

* ``"existing"``: start from the current problem state
* ``"zero"``: zero the scalar flux, precursor state, and stored angular fluxes
  before beginning the transient

.. note::

   ``initial_state="existing"`` is the usual choice when a transient begins
   from a previously computed steady-state or eigenvalue solution.

.. note::

   ``initial_state="zero"`` is mostly useful for controlled studies where the
   transient should begin from a clean zero field rather than from whatever
   state the problem currently holds.

Using Execute
-------------

The simplest transient pattern is:

.. code-block:: python

   phys.SetTimeDependentMode()

   solver = TransientSolver(
       problem=phys,
       dt=1.0e-3,
       stop_time=0.1,
       theta=BackwardEuler,
   )
   solver.Initialize()
   solver.Execute()

During :py:meth:`Execute`, the solver repeatedly advances until the current time
reaches ``stop_time``.

Important behavior:

* the last timestep is shortened automatically if needed so the solver lands
  exactly on ``stop_time``
* pre- and post-advance callbacks, if registered, are called around each
  internal timestep

.. note::

   This transient path assumes the problem was created with
   ``options={"save_angular_flux": True}``. Without that option, transient
   setup is invalid.

Use:

* :py:meth:`SetPreAdvanceCallback`
* :py:meth:`SetPostAdvanceCallback`

for work such as adaptive timestep changes, output, or diagnostics.

Example:

.. code-block:: python

   solver = TransientSolver(problem=phys, dt=1.0e-3, stop_time=0.05)
   solver.Initialize()
   ff = phys.GetScalarFluxFieldFunction()

   def post_advance():
       print("time =", phys.GetTime())
       for field in ff:
           field.Update()

   solver.SetPostAdvanceCallback(post_advance)
   solver.Execute()

.. note::

   The callback-based ``Execute`` path is the most convenient when the timestep
   loop is conceptually owned by the solver and Python is only observing or
   lightly adjusting the run.

.. note::

   A good rule is: if you can describe the transient as "run to this final time
   and do a little work each step," use ``Execute``. If you instead think "each
   next step depends on Python decisions made after the previous one," use a
   manual loop.

Using an Explicit Python Timestep Loop
--------------------------------------

You can also drive the transient manually with :py:meth:`Advance`.

This is useful when:

* the timestep size depends on Python-side logic,
* material maps or sources are being changed between steps,
* output logic is easier to express directly in Python, or
* you want explicit control over exactly when each timestep is taken

Pattern:

.. code-block:: python

   phys.SetTimeDependentMode()

   solver = TransientSolver(problem=phys, dt=1.0e-3, theta=CrankNicolson)
   solver.Initialize()
   ff = phys.GetScalarFluxFieldFunction()

   for step in range(num_steps):
       if step < 10:
           solver.SetTimeStep(1.0e-4)
       else:
           solver.SetTimeStep(5.0e-4)

       solver.Advance()

       print("time =", phys.GetTime())
       print("fission production =", phys.ComputeFissionProduction("new"))

       for field in ff:
           field.Update()

Important points:

* :py:meth:`Initialize` must be called before the first :py:meth:`Advance`
* :py:meth:`SetTimeStep` and :py:meth:`SetTheta` act on subsequent calls to
  :py:meth:`Advance`
* in a manual loop, *you* decide when to stop
* in a manual loop, ``stop_time`` is not the main control variable unless you
  choose to use it in your own Python logic
* transient problems still require ``options={"save_angular_flux": True}`` on
  the problem object

.. note::

   The explicit Python loop is usually the better design when the transient is
   being coupled to user logic, control rules, or parameter changes between
   timesteps. It keeps the sequencing visible in the script.

.. note::

   The two transient styles are complementary. ``Execute`` is better when the
   time march is conceptually fixed and Python only needs hooks. A manual
   ``Advance`` loop is better when Python is actively steering the run.

Balance Reporting
-----------------

:py:meth:`TransientSolver.ComputeBalanceTable` returns a global timestep balance
summary that includes:

* absorption, production, inflow, and outflow rates
* initial and final inventory
* predicted and actual inventory change
* inventory residual

This is often the most direct way to check that a transient step is behaving
sensibly.

.. note::

   When a transient result looks suspicious, the balance table is often the
   fastest sanity check because it summarizes whether the step is behaving like
   a physically consistent inventory update.

Power Iteration k-Eigen Solver
==============================

:py:class:`pyopensn.solver.PowerIterationKEigenSolver` solves the k-eigenvalue
problem with outer power iteration.

Use it when:

* you need a standard k-eigenvalue solve, and
* the usual power-iteration convergence behavior is acceptable

Constructor
-----------

The constructor takes:

* ``problem``: an existing
  :py:class:`pyopensn.solver.DiscreteOrdinatesProblem`
* ``max_iters``: maximum power iterations
* ``k_tol``: convergence tolerance on ``k_eff``
* ``reset_phi0``: whether to reset scalar fluxes to 1.0 before solving

GPU support
-----------

This solver can run with GPU acceleration when the associated
:py:class:`pyopensn.solver.DiscreteOrdinatesProblem` is configured for GPU
sweeps. The same problem-side restrictions apply as for steady-state source
solves: GPU use requires ``use_gpus=True`` with ``sweep_type="AAH"`` on a
supported non-curvilinear problem.

Example:

.. code-block:: python

   solver = PowerIterationKEigenSolver(
       problem=phys,
       max_iters=500,
       k_tol=1.0e-8,
   )
   solver.Initialize()
   solver.Execute()
   print(solver.GetEigenvalue())

How It Operates
---------------

At a high level, the solver:

1. optionally resets the starting scalar flux
2. forms the fission source from the current iterate
3. solves the transport problem with the problem's AGS/groupset iteration
4. updates ``k_eff``
5. repeats until the eigenvalue change satisfies ``k_tol`` or ``max_iters`` is
   reached

This solver uses the groupset inner solver and AGS settings configured on the
problem.

.. note::

   Power iteration is usually the first k-eigen method to try because it uses
   the same problem configuration style as steady-state source solves and is
   easy to reason about.

.. note::

   If a power-iteration run converges too slowly, first look at the transport
   iteration quality on the problem and the quality of the starting state before
   reaching for a different eigen solver.

.. note::

   In many workflows, power iteration is also used as the "reference" eigen
   solve because its iteration history is easy to interpret and compare against
   other methods.

Balance reporting with :py:meth:`ComputeBalanceTable` uses
``1 / k_eff`` scaling on the production term.

Nonlinear k-Eigen Solver
========================

:py:class:`pyopensn.solver.NonLinearKEigenSolver` solves the k-eigenvalue
problem with a nonlinear solve strategy rather than outer power iteration.

Use it when:

* you want a nonlinear k-eigen solve,
* you need its convergence behavior, or
* you want a small number of initial power iterations before entering the
  nonlinear solve

Constructor
-----------

The constructor takes:

* ``problem``: an existing :py:class:`pyopensn.solver.DiscreteOrdinatesProblem`
* nonlinear tolerances:
  ``nl_abs_tol``, ``nl_rel_tol``, ``nl_sol_tol``, ``nl_max_its``
* linear tolerances:
  ``l_abs_tol``, ``l_rel_tol``, ``l_div_tol``, ``l_max_its``
* GMRES controls:
  ``l_gmres_restart_intvl``, ``l_gmres_breakdown_tol``
* ``reset_phi0``
* ``num_initial_power_iterations``

GPU support
-----------

This solver does not currently provide a supported GPU execution path.

Example:

.. code-block:: python

   solver = NonLinearKEigenSolver(
       problem=phys,
       nl_abs_tol=1.0e-8,
       nl_max_its=50,
       l_abs_tol=1.0e-8,
       l_max_its=50,
       num_initial_power_iterations=3,
   )
   solver.Initialize()
   solver.Execute()
   print(solver.GetEigenvalue())

Important distinction:

* this solver does **not** use the groupset ``inner_linear_method`` setting for
  its internal nonlinear linearizations
* instead, it uses its own solver-level nonlinear and linear tolerance
  parameters
* the Python API exposes GMRES-specific controls for those internal linear
  solves, but does not expose an alternative internal linear method choice

.. note::

   This distinction is easy to miss. If a nonlinear k-eigen run needs tuning,
   look first at the ``nl_*`` and ``l_*`` options on
   :py:class:`pyopensn.solver.NonLinearKEigenSolver`, not at the groupset inner
   solver method.

.. warning::

   The nonlinear solver's internal GMRES controls have the same memory tradeoff
   as groupset GMRES. Larger ``l_gmres_restart_intvl`` values retain more Krylov
   vectors and can materially increase memory usage on large problems.

.. note::

   A practical way to choose between the two k-eigen solvers is to start with
   power iteration unless you already know you need the nonlinear method's
   behavior or controls. That keeps the initial setup simpler.

As with power iteration, :py:meth:`ComputeBalanceTable` uses
``1 / k_eff`` scaling on the production term.

Initialization Order and Common Patterns
========================================

A few patterns are worth making explicit.

Steady-state source solve:

.. code-block:: python

   phys = DiscreteOrdinatesProblem(...)
   solver = SteadyStateSourceSolver(problem=phys)
   solver.Initialize()
   solver.Execute()

Steady-state k-eigen solve:

.. code-block:: python

   phys = DiscreteOrdinatesProblem(...)
   solver = PowerIterationKEigenSolver(problem=phys)
   solver.Initialize()
   solver.Execute()

Steady-state to transient handoff:

.. code-block:: python

   phys = DiscreteOrdinatesProblem(...)

   ss = SteadyStateSourceSolver(problem=phys)
   ss.Initialize()
   ss.Execute()

   phys.SetTimeDependentMode()

   tr = TransientSolver(problem=phys, dt=1.0e-3, initial_state="existing")
   tr.Initialize()
   tr.Execute()

Manual transient loop:

.. code-block:: python

   phys.SetTimeDependentMode()
   tr = TransientSolver(problem=phys, dt=1.0e-3)
   tr.Initialize()

   while phys.GetTime() < 0.1:
       tr.Advance()

.. note::

   The most common initialization mistake is trying to use a solver against a
   problem that is in the wrong mode. Steady-state and k-eigen solvers require a
   steady-state problem; the transient solver requires a time-dependent problem.

.. note::

   The second most common mistake is forgetting that ``Initialize()`` is part of
   the contract. In OpenSn, initialization is not just a formality; it sets up
   solver state, restart behavior, and mode-dependent internals before the first
   solve step.

Field Functions and Solver Output
=================================

Problem objects, not solver objects, own the field-function accessors.

Common post-solve access patterns are:

* :py:meth:`pyopensn.solver.LBSProblem.GetScalarFluxFieldFunction`
* :py:meth:`pyopensn.solver.LBSProblem.CreateFieldFunction`
* :py:meth:`pyopensn.solver.DiscreteOrdinatesProblem.GetAngularFieldFunctionList`

For transient problems:

* if you are using :py:meth:`Execute`, query or export field functions in a
  post-advance callback
* if you are using a manual Python loop, query or export them after each call to
  :py:meth:`Advance`
* if you reuse field-function objects across solves or timesteps, call
  ``Update()`` before exporting or interpolating them

.. note::

   This is another reason to keep the distinction between problem and solver
   clear: the solver advances the state, but the problem is where you inspect
   the resulting transport fields.

.. note::

   In the API, these field functions are part of the problem-side output
   interface. That makes them easy to use from any solver, but it also means
   the decision to store certain outputs, such as angular flux, is made in the
   problem configuration.
