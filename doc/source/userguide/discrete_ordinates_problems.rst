===========================
Discrete-Ordinates Problems
===========================

The discrete-ordinates problem object is the center of OpenSn inputs. It gathers
the mesh, materials, groupsets, sources, boundary conditions, problem options,
and transport state into one object that can then be driven by a solver.

The two main problem classes exposed in Python are:

* :py:class:`pyopensn.solver.DiscreteOrdinatesProblem`
* :py:class:`pyopensn.solver.DiscreteOrdinatesCurvilinearProblem`

Overview
========

Use :py:class:`pyopensn.solver.DiscreteOrdinatesProblem` for Cartesian
discrete-ordinates problems. This is the standard problem class for most
steady-state, transient, and k-eigenvalue workflows.

The usual Cartesian construction pattern is:

.. code-block:: python

   phys = DiscreteOrdinatesProblem(
       mesh=mesh,
       num_groups=num_groups,
       groupsets=groupsets,
       xs_map=xs_map,
       boundary_conditions=boundary_conditions,
       point_sources=point_sources,
       volumetric_sources=volumetric_sources,
       options={...},
       sweep_type="AAH",
       use_gpus=False,
       time_dependent=False,
   )

The problem object owns:

* the mesh and material mapping,
* the scalar and angular transport state,
* the groupset definitions,
* the current forward or adjoint mode,
* the current steady-state or time-dependent mode,
* the field-function interface associated with the problem.

In practice, that means the problem is where the transport model is assembled:

* the spatial domain and geometry ids come from the mesh,
* the material assignment comes from ``xs_map``,
* the energy groupings and inner iteration setup come from ``groupsets``,
* driving terms come from sources and boundary conditions,
* global transport behavior comes from problem-level options.

.. note::

   A solver does not define the transport model. The problem does. The solver
   only decides how that already-defined model is advanced or converged.

Constructor Summary
===================

The main constructor inputs are:

* ``mesh``: a :py:class:`pyopensn.mesh.MeshContinuum`
* ``num_groups``: total number of energy groups
* ``groupsets``: one or more groupset dictionaries
* ``xs_map``: block-id to cross-section mapping
* ``boundary_conditions``: optional boundary specifications
* ``point_sources``: optional point sources
* ``volumetric_sources``: optional volumetric sources
* ``options``: problem-level settings
* ``sweep_type``: ``"AAH"`` or ``"CBC"``
* ``use_gpus``: request GPU sweep support where available
* ``time_dependent``: whether the problem starts in time-dependent mode

A complete problem definition often looks like:

.. code-block:: python

   phys = DiscreteOrdinatesProblem(
       mesh=mesh,
       num_groups=2,
       groupsets=groupsets,
       xs_map=[
           {"block_ids": [0], "xs": fuel_xs},
           {"block_ids": [1], "xs": moderator_xs},
       ],
       boundary_conditions=boundary_conditions,
       volumetric_sources=volumetric_sources,
       options={
           "verbose_inner_iterations": False,
           "verbose_outer_iterations": True,
       },
       sweep_type="AAH",
       time_dependent=False,
   )

Each of these pieces is discussed elsewhere in the guide, but it is useful to
see them in one place.

Constructor Inputs
==================

``mesh``
--------

This is the transport mesh.

It must already exist and already carry the block ids and boundary labels that
the rest of the problem will use.

``num_groups``
--------------

This is the total number of energy groups in the problem.

It must be consistent with:

* the loaded cross sections,
* the groupset definitions,
* any group-wise source data.

``groupsets``
-------------

This is the list of groupset dictionaries that define angular quadratures,
inner iteration, angle aggregation, and optional DSA configuration.

Groupsets are detailed separately in :doc:`groupsets`, but from the problem
perspective they define how the energy range is partitioned for transport
iterations.

``xs_map``
----------

This maps mesh block ids to cross-section objects.

Every block id used in the transport region should be assigned an XS object.

``boundary_conditions``
-----------------------

This is the list of transport boundary-condition dictionaries. The names used
here must match the boundary ids on the mesh.

``point_sources`` and ``volumetric_sources``
--------------------------------------------

These are optional source objects added at construction time.

``options``
-----------

This dictionary holds problem-level behavior that is not naturally part of the
mesh, source, or groupset definitions.

``sweep_type``
--------------

This selects the sweep type:

* ``"AAH"``: the aggregated sweeper with cycle breaking
* ``"CBC"``: the cell-by-cell sweeper

``use_gpus``
------------

This requests GPU acceleration where supported.

``time_dependent``
------------------

This puts the problem into time-dependent mode at construction time.

Shared LBS Interface
====================

Because :py:class:`DiscreteOrdinatesProblem` derives from
:py:class:`pyopensn.solver.LBSProblem`, it also exposes the common LBS problem
methods such as:

* :py:meth:`GetScalarFluxFieldFunction`
* :py:meth:`CreateFieldFunction`
* :py:meth:`ComputeFissionRate`
* :py:meth:`WriteFluxMoments`
* :py:meth:`ReadFluxMoments`
* :py:meth:`SetPointSources`
* :py:meth:`SetVolumetricSources`
* :py:meth:`SetBoundaryOptions`
* :py:meth:`SetAdjoint`

Problem Options
===============

The Python API exposes problem options through the constructor and
:py:meth:`SetOptions`.

The most important options for everyday use are:

* ``adjoint``
* ``save_angular_flux``
* ``verbose_inner_iterations``
* ``verbose_outer_iterations``
* ``max_ags_iterations``
* ``ags_tolerance``
* ``ags_convergence_check``
* ``field_function_prefix_option``
* ``field_function_prefix``

There are also restart, precursor, and message-size options for more
specialized workflows.

.. note::

   Problem options are where users should look for global transport behavior.
   If a setting changes how the entire problem behaves rather than how one
   groupset behaves, it usually belongs here.

Setting options after construction
----------------------------------

Problem options can also be changed after construction:

.. code-block:: python

   phys.SetOptions(
       verbose_inner_iterations=False,
       max_ags_iterations=200,
       ags_tolerance=1.0e-8,
   )

This is useful for parameter studies or staged workflows where the problem
definition is reused.

``sweep_type``
==============

The Python API exposes:

* ``"AAH"``
* ``"CBC"``

Example:

.. code-block:: python

   sweep_type="AAH"

If omitted, the default is ``"AAH"``.

Operationally, the two sweep types are different:

* ``"AAH"`` is the more general aggregated sweeper and should be treated as the
  default production choice.
* ``"CBC"`` is a cell-by-cell sweeper that preserves exact cell-to-cell
  dependencies.

The practical differences are:

* ``AAH`` has explicit delayed-angular-flux machinery for cycle handling.
* ``AAH`` can break both local and inter-partition sweep cycles by removing
  feedback-arc-set edges and lagging the corresponding angular-flux
  dependencies.
* ``CBC`` does not support local sweep cycles.

.. note::

   In the ``AAH`` implementation, the lagged data is tied specifically to
   cycle-breaking dependencies. It is not a general statement that all angular
   fluxes are always lagged.

Practical recommendation:

* ``AAH`` remains the default production choice and is the safer option for
  most users, particularly for problems with cyclic sweep dependencies.
* Both ``AAH`` and ``CBC`` support time-dependent (transient) mode.
* Choose ``CBC`` only when the sweep graph is known to be acyclic or when
  you have verified it meets the acyclicity requirement for your specific
  problem.

Problem Modes
=============

``time_dependent``
------------------

If ``time_dependent=True``, the problem starts in time-dependent mode.

Important requirement:

* time-dependent operation requires ``options={"save_angular_flux": True}``

Example:

.. code-block:: python

   phys = DiscreteOrdinatesProblem(
       ...,
       time_dependent=True,
       options={"save_angular_flux": True},
   )

This requirement exists because transient updates need access to angular-flux
state from one timestep to the next.

This option changes the problem mode, not just a solver setting. That is why
it belongs on the problem object and must be consistent with the solver used
later.

``use_gpus``
------------

``use_gpus`` requests GPU acceleration for supported sweep paths.

Current restrictions:

* only ``"AAH"`` is supported for GPU use,
* curvilinear problems do not support GPU acceleration,
* time-dependent problems do not support GPU acceleration,
* adjoint problems do not support GPU acceleration.

Most users should treat this as a deployment choice after the base problem is
already running correctly on the CPU.

Adjoint Mode
------------

Adjoint mode is controlled at the problem level because it changes the meaning
of the transport problem itself, not just the iterative algorithm.

It can be set at construction time through options or later with
``SetAdjoint(True)``.

Because this is a fundamental change in problem interpretation, users should be
deliberate about when they switch it on.

Field-Function Interface
========================

The problem object is also where users access transport outputs.

The main field-function accessors are:

* :py:meth:`GetScalarFluxFieldFunction`
* :py:meth:`CreateFieldFunction`
* :py:meth:`GetAngularFieldFunctionList`

Example:

.. code-block:: python

   scalar_ffs = phys.GetScalarFluxFieldFunction()
   power_ff = phys.CreateFieldFunction("power_generation", "power")

.. note::

   In LBS workflows, users normally work with the field-function objects
   returned by these accessors directly. They are created from the current
   transport state when requested. If the transport state changes later, call
   ``Update()`` on an existing updateable field function or create a fresh field
   function from the current state.

Angular-Flux Access
===================

The discrete-ordinates problem also exposes angular-flux data directly through
:py:meth:`GetPsi`.

``GetPsi()`` returns a list of NumPy arrays that view the underlying angular-
flux storage.

This is a specialized interface and is mainly useful when:

* a workflow needs direct access to angular-flux data,
* a custom analysis step is easier to write in Python,
* the problem was configured with ``save_angular_flux=True`` for the needed
  workflow.

.. note::

   Angular flux is much larger than scalar flux. Many workflows never need it.
   Only enable and use it when the calculation actually requires it.

Balance and Leakage
===================

The Cartesian discrete-ordinates problem exposes two important diagnostics:

* :py:meth:`ComputeBalance`
* :py:meth:`ComputeLeakage`

``ComputeBalance()``
--------------------

This computes the particle balance for the problem and returns a dictionary of
balance terms.

It is useful for:

* checking whether the run is physically consistent,
* confirming source, absorption, leakage, and production trends,
* validating test and regression problems.

``ComputeLeakage()``
--------------------

This computes boundary leakage and returns a dictionary mapping boundary names
to group-wise NumPy arrays.

Example:

.. code-block:: python

   leakage = phys.ComputeLeakage(["xmin", "xmax"])

Important requirement:

* leakage computation requires ``save_angular_flux=True``

This requirement exists because leakage is derived from the outgoing angular
flux.

Writing and Reading Transport State
===================================

The problem object also exposes file-based state helpers such as:

* :py:meth:`WriteFluxMoments`
* :py:meth:`ReadFluxMoments`
* :py:meth:`CreateAndWriteSourceMoments`
* :py:meth:`ReadSourceMoments`
* :py:meth:`ReadFluxMomentsAndMakeSourceMoments`
* :py:meth:`WriteAngularFluxes`
* :py:meth:`ReadAngularFluxes`

These are useful for:

* restart-like workflows,
* response studies,
* source-driven workflows that reuse previously computed fields.

Updating the Problem In Place
=============================

Several parts of the problem can be updated after construction.

The most important methods are:

* :py:meth:`SetOptions`
* :py:meth:`SetPointSources`
* :py:meth:`SetVolumetricSources`
* :py:meth:`SetBoundaryOptions`
* :py:meth:`SetAdjoint`
* :py:meth:`SetTimeDependentMode`

Updating a Problem After Construction
-------------------------------------

One of the strengths of the problem API is that major parts of the model can be
replaced without rebuilding the entire problem object.

Examples:

.. code-block:: python

   phys.SetVolumetricSources(
       clear_volumetric_sources=True,
       volumetric_sources=[new_source],
   )

   phys.SetBoundaryOptions(
       clear_boundary_conditions=True,
       boundary_conditions=[new_boundary],
   )

   phys.SetAdjoint(True)

This makes the problem object usable in:

* source studies,
* boundary-condition studies,
* forward/adjoint comparisons,
* transient driver loops that change forcing terms.

DiscreteOrdinatesCurvilinearProblem
===================================

:py:class:`pyopensn.solver.DiscreteOrdinatesCurvilinearProblem` is the
curvilinear companion to the Cartesian problem class.

It uses the same general construction pattern, but currently requires:

* a suitable curvilinear mesh,
* ``coord_system=2`` for cylindrical coordinates,
* a compatible quadrature and solver setup.

Important current limitations:

* the curvilinear solver only supports cylindrical geometries,
* GPU acceleration is not supported,
* users should treat it as a more specialized path than the standard Cartesian
  problem.

Example:

.. code-block:: python

   phys = DiscreteOrdinatesCurvilinearProblem(
       mesh=mesh,
       coord_system=2,
       num_groups=num_groups,
       groupsets=groupsets,
       xs_map=xs_map,
       sweep_type="AAH",
   )

Typical Construction Patterns
=============================

A few patterns show up repeatedly:

* source-driven steady-state problem plus ``SteadyStateSourceSolver``
* time-dependent problem plus ``TransientSolver``
* multiplication problem plus ``PowerIterationKEigenSolver``
* multiplication problem plus ``NonLinearKEigenSolver``

Practical Guidance
==================

As a rule:

* use :py:class:`DiscreteOrdinatesProblem` unless the curvilinear class is
  specifically needed,
* keep ``sweep_type="AAH"`` unless there is a clear reason to choose ``CBC``,
* treat ``save_angular_flux`` as a capability switch with a memory cost,
* build the simplest correct problem first, then add extra options,
* remember that the problem object remains central even after the solver is
  constructed because it owns the transport state and output interface.

.. note::

   Many solver issues are actually problem-definition issues. If something looks
   unstable or physically wrong, it is often worth reviewing the problem object
   first: mesh labels, ``xs_map``, groupsets, source definitions, and problem
   options.
