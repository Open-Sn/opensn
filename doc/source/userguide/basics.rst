=============
OpenSn Basics
=============

This section is the short version of the user guide. If you only remember a few
things about writing OpenSn inputs, remember these.

The Core Model
==============

An OpenSn transport input is usually built from six pieces:

1. a mesh,
2. block ids and boundary names on that mesh,
3. cross sections assigned to block ids,
4. one or more groupsets,
5. a problem object,
6. a solver object.

In other words:

* the mesh defines the geometry,
* ``xs_map`` defines the materials,
* ``groupsets`` define the engery group structure and iteration parameters,
* the problem object defines the physics model,
* the solver object defines how that model is advanced.

.. note::

   The problem owns the physics model. The solver owns the algorithm. Keeping
   those roles separate makes most of the rest of the API easier to understand.

The Typical Workflow
====================

Most Python inputs follow this pattern:

.. code-block:: python

   grid = ...
   xs = ...
   groupsets = ...

   phys = DiscreteOrdinatesProblem(
       mesh=grid,
       num_groups=...,
       groupsets=groupsets,
       xs_map=[...],
       boundary_conditions=[...],
       volumetric_sources=[...],
       options={...},
   )

   solver = SteadyStateSourceSolver(problem=phys)
   solver.Initialize()
   solver.Execute()

   scalar_ffs = phys.GetScalarFluxFieldFunction()

If the problem is transient, use :py:class:`TransientSolver`.
If it is a criticality problem, use one of the k-eigen solvers.

What Usually Goes Where
=======================

As a rule of thumb:

* boundary names and block ids belong on the mesh,
* materials belong in ``xs_map``,
* angular quadrature and inner iteration parameters belong in ``groupsets``,
* outer-iteration parameters belong in the solver.

If you are not sure where an option belongs, ask whether it changes:

* a groupset solve,
* the whole problem, or
* the outer solve.

That usually tells you the right place.

The Three Most Important Practical Rules
========================================

1. Start simple.

   Start with one groupset, simple boundaries, and a small mesh if possible.

2. Converge the transport solve properly.

   If you have multiple groupsets and there is across-groupset upscatter, enable
   and converge the across-groupset solver (AGS). Do not rely on loose inner
   solves to be corrected by outer iterations later.

3. Postprocess through field functions.

   The usual output path is to solve first, then get scalar or power field
   functions and export or interpolate them.

.. note::

   Many user problems are easier to diagnose if you first ask: is this a mesh
   issue, a material/block id issue, a groupset/iteration issue, or an outer
   solver issue?

Good Default Habits
===================

Good starting habits are:

* use one groupset unless you know why you need more,
* use ``petsc_gmres`` as the default difficult-problem inner method,
* use ``angle_aggregation_type="single"`` on unstructured meshes,
* keep inner tolerances tighter than outer tolerances,
* create power field functions only when you need them,
* inspect exported field functions or interpolated values early in a new model.

.. note::

   ``petsc_gmres`` is a strong default for difficult problems, but it stores
   Krylov basis vectors between restarts. On large problems, memory usage can
   become significant, so keep ``gmres_restart_interval`` only as large as the
   problem needs.

What To Read Next
=================

If you want the next level of detail:

* read :doc:`geometry_mesh` for mesh creation and labeling,
* read :doc:`materials_xs` for materials and cross sections,
* read :doc:`groupsets` and :doc:`iterative_methods` for transport iteration,
* read :doc:`discrete_ordinates_problems` and :doc:`solvers` for problem and
  solver setup,
* read :doc:`example_problems` for complete input patterns,
* read :doc:`post_processors` for field functions and interpolation.
