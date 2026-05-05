================
Example Problems
================

This section collects annotated example-input patterns for many of the
problem-and-solver combinations used in OpenSn. The examples are intentionally
verbose and repetitive. They are meant to show where each piece of an input
belongs rather than to be as short as possible.

The goal is that a user should be able to read this section and understand how
to assemble a complete transport input from mesh generation through field-
function output.

Common Input Structure
======================

Most OpenSn inputs follow the same high-level structure:

1. create or import the mesh,
2. assign block ids and boundary names,
3. create or load cross sections,
4. define the ``xs_map``,
5. define one or more groupsets,
6. define optional sources and boundary conditions,
7. construct the problem object,
8. construct the solver object,
9. initialize and execute,
10. retrieve or export field functions.

Overview
========

In practice:

* the mesh section defines the geometry and ids used later,
* the material section defines what data lives on which block ids,
* groupsets hold the angular quadrature and iterative controls,
* sources and boundaries define the forcing and transport closure,
* the problem object collects the physical model,
* the solver object selects the algorithm,
* field functions provide the usual output path.

.. note::

   A useful way to read the examples below is to separate them into
   "physics-definition" blocks and "algorithm" blocks. The problem object holds
   the physical model; the solver object decides how that model is advanced.

Common Building Blocks
======================

Before looking at the full examples, it helps to identify the pieces that
appear in almost every script.

Mesh
----

.. code-block:: python

   from pyopensn.mesh import OrthogonalMeshGenerator

   mesh = OrthogonalMeshGenerator(
       node_sets=[
           [0.0, 1.0, 2.0, 3.0],
           [0.0, 1.0],
       ],
   ).Execute()
   mesh.SetOrthogonalBoundaries()
   mesh.SetUniformBlockID(0)

Cross Sections
--------------

.. code-block:: python

   from pyopensn.xs import MultiGroupXS

   xs = MultiGroupXS()
   xs.LoadFromOpenSn("material_2g.cxs")

   xs_map = [{"block_ids": [0], "xs": xs}]

Groupsets
---------

.. code-block:: python

   from pyopensn.aquad import GLCProductQuadrature2DXY

   groupsets = [
       {
           "groups_from_to": [0, 1],
           "angular_quadrature": GLCProductQuadrature2DXY(
               n_azimuthal=8,
               n_polar=4,
               scattering_order=1,
           ),
           "inner_linear_method": "petsc_gmres",
           "l_abs_tol": 1.0e-6,
           "l_max_its": 200,
           "gmres_restart_interval": 30,
       }
   ]

Boundary Conditions
-------------------

.. code-block:: python

   boundary_conditions = [
       {"name": "xmin", "type": "vacuum"},
       {"name": "xmax", "type": "vacuum"},
       {"name": "ymin", "type": "reflecting"},
       {"name": "ymax", "type": "reflecting"},
   ]

Volumetric Sources
------------------

.. code-block:: python

   from pyopensn.source import VolumetricSource

   volumetric_sources = [
       VolumetricSource(
           block_ids=[0],
           group_strength=[1.0, 0.0],
       )
   ]

Point Sources
-------------

.. code-block:: python

   from pyopensn.source import PointSource

   point_sources = [
       PointSource(
           location=(0.5, 0.5, 0.0),
           strength=[1.0, 0.0],
       )
   ]

Problem Construction
--------------------

.. code-block:: python

   from pyopensn.solver import DiscreteOrdinatesProblem

   phys = DiscreteOrdinatesProblem(
       mesh=mesh,
       num_groups=2,
       groupsets=groupsets,
       xs_map=xs_map,
       boundary_conditions=boundary_conditions,
       volumetric_sources=volumetric_sources,
       options={
           "verbose_inner_iterations": False,
           "verbose_outer_iterations": True,
       },
   )

.. note::

   Most user scripts spend more effort building the problem than building the
   solver. That is normal. The problem object captures the physical model.

Example 1: Source-Driven Steady-State Problem
=============================================

This is the standard starting point for most users.

.. code-block:: python

   import pyopensn as opensn
   from pyopensn.mesh import OrthogonalMeshGenerator
   from pyopensn.xs import MultiGroupXS
   from pyopensn.aquad import GLCProductQuadrature2DXY
   from pyopensn.source import VolumetricSource
   from pyopensn.solver import DiscreteOrdinatesProblem, SteadyStateSourceSolver
   from pyopensn.fieldfunc import FieldFunctionGridBased

   mesh = OrthogonalMeshGenerator(
       node_sets=[
           [0.0, 1.0, 2.0, 3.0],
           [0.0, 1.0, 2.0],
       ],
   ).Execute()
   mesh.SetOrthogonalBoundaries()
   mesh.SetUniformBlockID(0)

   xs = MultiGroupXS()
   xs.LoadFromOpenSn("material_2g.cxs")

   quadrature = GLCProductQuadrature2DXY(
       n_azimuthal=8,
       n_polar=4,
       scattering_order=1,
   )

   groupsets = [
       {
           "groups_from_to": [0, 1],
           "angular_quadrature": quadrature,
           "inner_linear_method": "petsc_gmres",
           "l_abs_tol": 1.0e-6,
           "l_max_its": 200,
       }
   ]

   xs_map = [
       {"block_ids": [0], "xs": xs},
   ]

   boundary_conditions = [
       {"name": "xmin", "type": "vacuum"},
       {"name": "xmax", "type": "vacuum"},
       {"name": "ymin", "type": "vacuum"},
       {"name": "ymax", "type": "vacuum"},
   ]

   volumetric_sources = [
       VolumetricSource(block_ids=[0], group_strength=[1.0, 0.0]),
   ]

   phys = DiscreteOrdinatesProblem(
       mesh=mesh,
       num_groups=2,
       groupsets=groupsets,
       xs_map=xs_map,
       boundary_conditions=boundary_conditions,
       volumetric_sources=volumetric_sources,
       options={
           "verbose_inner_iterations": False,
           "verbose_outer_iterations": True,
       },
   )

   solver = SteadyStateSourceSolver(problem=phys)
   solver.Initialize()
   solver.Execute()

   scalar_ffs = phys.GetScalarFluxFieldFunction()
   FieldFunctionGridBased.ExportMultipleToPVTU(
       scalar_ffs,
       "steady_state_flux",
   )

What goes where:

* ``mesh`` contains the geometry and block/boundary labels,
* ``xs_map`` assigns materials to block ids,
* ``groupsets`` define quadrature, sweep aggregation, and inner iteration,
* ``boundary_conditions`` define the transport behavior at the mesh boundary,
* ``volumetric_sources`` and ``point_sources`` provide fixed driving terms,
* ``options`` contain problem-wide settings such as verbosity and AGS controls,
* the steady-state solver provides the outer solve algorithm.

.. note::

   When debugging a steady-state problem, it is usually best to start with a
   single groupset, low scattering order, and a simple vacuum or reflecting
   boundary setup. Once the baseline case behaves correctly, then add more
   detail.

Example 2: Multi-Material Fixed-Source Problem
==============================================

The next common step is assigning different cross sections to different mesh
regions.

.. code-block:: python

   fuel_xs = MultiGroupXS()
   fuel_xs.LoadFromOpenSn("fuel.cxs")

   moderator_xs = MultiGroupXS()
   moderator_xs.LoadFromOpenSn("moderator.cxs")

   mesh.SetBlockIDFromLogicalVolume(fuel_lv, 0, True)
   mesh.SetBlockIDFromLogicalVolume(moderator_lv, 1, True)

   xs_map = [
       {"block_ids": [0], "xs": fuel_xs},
       {"block_ids": [1], "xs": moderator_xs},
   ]

The rest of the problem definition can remain unchanged.

.. note::

   A large fraction of practical "geometry work" in an OpenSn input is really
   about assigning the right block ids to the right cells. That is why block id
   checking is worth doing early.

Example 3: Boundary-Driven Problem
==================================

Some problems are driven entirely or primarily by incoming boundary flux.

Isotropic inflow:

.. code-block:: python

   boundary_conditions = [
       {"name": "xmin", "type": "isotropic", "group_strength": [1.0, 0.0]},
       {"name": "xmax", "type": "vacuum"},
       {"name": "ymin", "type": "reflecting"},
       {"name": "ymax", "type": "reflecting"},
   ]

Arbitrary angular inflow:

.. code-block:: python

   from pyopensn.math import AngularFluxFunction

   def incident_flux(group, direction_num):
       if group == 0:
           return 1.0
       return 0.0

   boundary_conditions = [
       {
           "name": "xmin",
           "type": "arbitrary",
           "function": AngularFluxFunction(incident_flux),
       },
       {"name": "xmax", "type": "vacuum"},
   ]

.. note::

   Boundary-condition dictionaries only affect boundaries whose ids already
   exist on the mesh.

Example 4: Time-Dependent Source Problem
========================================

For a transient problem that does not require Python-side updates between
timesteps, ``Execute()`` is the simplest path.

.. code-block:: python

   from pyopensn.solver import TransientSolver, CrankNicolson

   phys = DiscreteOrdinatesProblem(
       mesh=mesh,
       num_groups=1,
       groupsets=groupsets,
       xs_map=xs_map,
       volumetric_sources=volumetric_sources,
       options={
           "save_angular_flux": True,
           "verbose_inner_iterations": False,
       },
       time_dependent=True,
   )

   solver = TransientSolver(
       problem=phys,
       dt=1.0e-3,
       stop_time=1.0e-2,
   )
   solver.Initialize()
   solver.SetTheta(CrankNicolson)
   solver.Execute()

Important changes relative to a steady-state solve:

* the problem is created with ``time_dependent=True``,
* ``save_angular_flux=True``,
* the outer solver becomes ``TransientSolver``,
* timestep-related settings belong to the solver, not the problem.

.. important::

   For transient problems, ``save_angular_flux=True`` is required. It is not
   just an optional post-processing switch.

Example 5: Explicit Python Time Loop
====================================

If the input needs to update the problem between timesteps, use ``Advance()``
directly.

.. code-block:: python

   from pyopensn.solver import TransientSolver, BackwardEuler

   dt = 1.0e-3
   stop_time = 1.0e-2
   t = 0.0

   solver = TransientSolver(
       problem=phys,
       dt=dt,
       stop_time=stop_time,
   )
   solver.Initialize()
   solver.SetTheta(BackwardEuler)
   scalar_ff = phys.GetScalarFluxFieldFunction()[0]

   while t < stop_time:
       solver.Advance()
       t += dt

       scalar_ff.Update()

       if t > 5.0e-3:
           dt = 5.0e-4
           solver.SetTimeStep(dt)

This pattern is useful when:

* timestep size changes during the run,
* a custom output cadence is needed,
* source or boundary data are updated from Python between steps.

.. note::

   The same transient requirement still applies in a manual loop:
   ``save_angular_flux=True`` must already be enabled on the problem before the
   transient solve begins.

.. note::

   Field functions are created from the current completed state when requested.
   If the state changes on a later step, call ``Update()`` on an updateable
   field function or create a fresh field function from the current state.

Example 6: Power-Iteration k-Eigenvalue Problem
===============================================

Power iteration is the usual first choice for a standard k-eigenvalue solve.

.. code-block:: python

   from pyopensn.solver import PowerIterationKEigenSolver

   phys = DiscreteOrdinatesProblem(
       mesh=mesh,
       num_groups=1,
       groupsets=groupsets,
       xs_map=xs_map,
       options={
           "power_default_kappa": 1.0,
       },
   )

   solver = PowerIterationKEigenSolver(
       problem=phys,
       max_iters=1000,
       k_tol=1.0e-10,
       reset_solution=True,
       reset_phi0=True,
   )
   solver.Initialize()
   solver.Execute()

   scalar_ffs = phys.GetScalarFluxFieldFunction()
   power_ff = phys.CreateFieldFunction(
       "power_generation_norm",
       "power",
       power_normalization_target=1.0,
   )

This is the standard pattern for k-eigenvalue work when power iteration is the
desired outer method and a power field is wanted for post processing.

Example 7: Nonlinear k-Eigenvalue Problem
=========================================

The nonlinear solver is configured through its own solver-level tolerances.

.. code-block:: python

   from pyopensn.solver import NonLinearKEigenSolver

   solver = NonLinearKEigenSolver(
       problem=phys,
       nl_abs_tol=1.0e-8,
       nl_rel_tol=1.0e-8,
       nl_sol_tol=1.0e-50,
       nl_max_its=50,
       l_abs_tol=1.0e-8,
       l_rel_tol=1.0e-8,
       l_div_tol=1.0e6,
       l_max_its=50,
       l_gmres_restart_intvl=30,
       l_gmres_breakdown_tol=1.0e6,
       reset_phi0=True,
       num_initial_power_iterations=2,
   )
   solver.Initialize()
   solver.Execute()

.. note::

   The nonlinear k-eigen solver has its own linear solver controls. These are
   separate from the groupset ``inner_linear_method`` used in the standard
   transport solve path. Groupset linear solver options are not used with the
   nonlinear k-eigen solver.

Example 8: Curvilinear Problem
==============================

Curvilinear problems follow the same general pattern, but use the curvilinear
problem type and a compatible quadrature.

.. code-block:: python

   from pyopensn.aquad import GLCProductQuadrature2DRZ
   from pyopensn.solver import (
       DiscreteOrdinatesCurvilinearProblem,
       SteadyStateSourceSolver,
   )

   quadrature = GLCProductQuadrature2DRZ(
       n_azimuthal=16,
       n_polar=8,
       scattering_order=1,
   )

   phys = DiscreteOrdinatesCurvilinearProblem(
       mesh=mesh,
       coord_system=2,
       num_groups=2,
       groupsets=[
           {
               "groups_from_to": [0, 1],
               "angular_quadrature": quadrature,
           }
       ],
       xs_map=xs_map,
   )

   solver = SteadyStateSourceSolver(problem=phys)
   solver.Initialize()
   solver.Execute()

The curvilinear problem supports CPU ``AAH`` and ``CBC`` sweep types. Use
``AAH`` as the default, particularly when cyclic sweep dependencies are
possible; choose ``CBC`` only when the sweep graph satisfies CBC's acyclicity
requirements.

Example 9: Updating a Problem In Place
======================================

Problem objects are mutable. An input can replace sources, boundaries, or
adjoint mode without rebuilding the entire problem object.

.. code-block:: python

   from pyopensn.source import VolumetricSource

   phys.SetVolumetricSources(
       clear_volumetric_sources=True,
       volumetric_sources=[
           VolumetricSource(
               block_ids=[0],
               group_strength=[0.1, 0.0],
           )
       ],
   )

   phys.SetBoundaryOptions(
       clear_boundary_conditions=True,
       boundary_conditions=[
           {"name": "xmin", "type": "vacuum"},
           {"name": "xmax", "type": "reflecting"},
       ],
   )

   phys.SetAdjoint(True)

This kind of in-place update is useful for source studies, transient driver
loops, and forward/adjoint comparison workflows.

Example 10: Field-Function Output
=================================

Post processing is often only a few lines, but it belongs at the end of the
workflow after a successful solve.

.. code-block:: python

   from pyopensn.fieldfunc import FieldFunctionGridBased

   scalar_ffs = phys.GetScalarFluxFieldFunction()
   FieldFunctionGridBased.ExportMultipleToPVTU(
       scalar_ffs,
       "phi",
   )

   if want_power:
       power_ff = phys.CreateFieldFunction("power_generation", "power")
       FieldFunctionGridBased.ExportMultipleToPVTU(
           [power_ff],
           "power",
       )

You can also export both together:

.. code-block:: python

   FieldFunctionGridBased.ExportMultipleToPVTU(
       [scalar_ffs[0], power_ff],
       "phi_and_power",
   )

Example 11: Balance and Leakage Checks
======================================

Post-solve diagnostics are often as important as the field export itself.

.. code-block:: python

   phys.ComputeBalance()

   leakage = phys.ComputeLeakage(["xmin", "xmax"])
   print(leakage["xmin"])

This is especially useful in regression tests and in the early validation of a
new physical model.

Example 12: Reusing Flux Data
=============================

Some workflows use saved flux moments or angular fluxes as inputs to later
stages:

.. code-block:: python

   phys.WriteFluxMoments("phi_dump")
   phys.ReadFluxMoments("phi_dump", single_file_flag=True)

or, when angular-flux storage is enabled:

.. code-block:: python

   phys.WriteAngularFluxes("psi_dump")
   phys.ReadAngularFluxes("psi_dump")

Example 13: Writing and Reading Restart Dumps
=============================================

Restart dumps are different from regular field-function or angular i/o. They
save a restartable state owned primarily by the problem, with small
solver-specific data added only when needed:

* common problem-owned state includes flux moments, time-integration metadata,
  precursor state, and any discrete-ordinates angular restart data required by
  the active problem mode,
* :py:class:`SteadyStateSourceSolver` uses only the common problem-owned
  restart state,
* :py:class:`PowerIterationKEigenSolver` adds k-eigenvalue-specific quantities
  such as ``k_eff`` and the previous fission-production normalization,
* :py:class:`TransientSolver` uses the common problem-owned restart state and
  adds the transient step counter.

To write restart dumps for a steady-state solve:

.. code-block:: python

   phys = DiscreteOrdinatesProblem(
       mesh=mesh,
       num_groups=num_groups,
       groupsets=groupsets,
       xs_map=xs_map,
       options={
           "restart_writes_enabled": True,
           "write_restart_path": "restart_data/my_problem",
           "write_delayed_psi_to_restart": True,
       },
   )

   solver = SteadyStateSourceSolver(problem=phys)
   solver.Initialize()
   solver.Execute()

Important points:

* ``restart_writes_enabled=True`` enables restart output,
* ``write_restart_path`` is the file stem to write,
* ``write_delayed_psi_to_restart=True`` includes delayed angular-flux state in
  the restart dump.

.. note::

   For :py:class:`SteadyStateSourceSolver`, the restart dump is written after
   the transport solve has completed. It is not written mid-iteration.

For long-running problems, timed restart writes can also be enabled:

.. code-block:: python

   options={
       "restart_writes_enabled": True,
       "write_restart_path": "restart_data/my_problem",
       "write_restart_time_interval": 60,
   }

This timed path is relevant to solver types that execute over multiple outer
iterations or time steps, including:
:py:class:`pyopensn.solver.PowerIterationKEigenSolver` and
:py:class:`pyopensn.solver.TransientSolver`. They check the elapsed wall-clock
time against ``write_restart_time_interval`` during execution and still write a
final restart dump at the end when restart output is enabled.

Then a later run can read that restart state:

.. code-block:: python

   phys = DiscreteOrdinatesProblem(
       mesh=mesh,
       num_groups=num_groups,
       groupsets=groupsets,
       xs_map=xs_map,
       options={
           "read_restart_path": "restart_data/my_problem",
       },
   )

   solver = SteadyStateSourceSolver(problem=phys)
   solver.Initialize()
   solver.Execute()

.. note::

   Restart reads happen during :py:meth:`Initialize`, not during
   :py:meth:`Execute`.

This is useful when:

* a large steady-state solve needs to be resumed,
* a k-eigen solve should restart from a previous state,
* a transient solve should continue from a previously written time step,
* or a workflow is split across multiple runs.

.. note::

   Restart reads and writes assume a compatible problem definition. In
   practice, that means the mesh, group structure, discretization, and related
   state layout must match the restart data being read.

Choosing a Template
===================

As a practical guide:

* start from Example 1 for a basic fixed-source calculation,
* use Example 3 if the transient has a fixed timestep structure,
* use Example 4 if Python must intervene between timesteps,
* use Example 5 for standard k-eigenvalue work,
* use Example 6 only when the nonlinear eigenvalue solve is specifically
  needed,
* use Example 13 when a workflow needs restartable transport state across runs,
* add the post-processing patterns only after the solve itself is behaving
  correctly.
