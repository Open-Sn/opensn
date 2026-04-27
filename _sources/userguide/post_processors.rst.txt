===============
Post Processors
===============

For most OpenSn workflows, post processing starts with field functions. A field
function is a mesh-based representation of transport data that can be exported
or interpolated from Python.

The important practical point is that field functions are created from the
current solver state when requested. They are snapshots, not continuously
updated solver-owned views. Some field functions can refresh that same
snapshot object explicitly through
:py:meth:`pyopensn.fieldfunc.FieldFunctionGridBased.Update`.

Overview
========

The OpenSn field function interfaces are:

* :py:meth:`pyopensn.solver.LBSProblem.GetScalarFluxFieldFunction`
* :py:meth:`pyopensn.solver.LBSProblem.CreateFieldFunction`
* :py:meth:`pyopensn.solver.DiscreteOrdinatesProblem.GetAngularFieldFunctionList`
* :py:meth:`pyopensn.fieldfunc.FieldFunctionGridBased.ExportMultipleToPVTU`
* :py:class:`pyopensn.fieldfunc.FieldFunctionInterpolationPoint`
* :py:class:`pyopensn.fieldfunc.FieldFunctionInterpolationLine`
* :py:class:`pyopensn.fieldfunc.FieldFunctionInterpolationVolume`

In practice, most workflows are:

1. solve the problem,
2. create the field functions needed from the current state,
3. export them or evaluate them with interpolation objects.

.. note::

   If the solver state changes, a field function created earlier does not
   update itself automatically after a later solve or timestep. Either create a
   new field function from the current state, or call ``Update()`` on the
   existing field function when ``CanUpdate()`` returns ``True``.

.. note::

   For LBS workflows, the important field-function creators are
   :py:meth:`GetScalarFluxFieldFunction`,
   :py:meth:`CreateFieldFunction`, and, for discrete ordinates problems,
   :py:meth:`GetAngularFieldFunctionList`. These methods return fresh field
   functions from the current state; they do not rely on a persistent
   field-function cache.

Updating Existing Field Functions
=================================

Field functions returned by LBS problem accessors are updateable snapshots.
They hold their own field-vector data, but also know how to refresh that data
from the owning problem while the owning problem is still alive.

Use :py:meth:`pyopensn.fieldfunc.FieldFunctionGridBased.CanUpdate` before
refreshing a field function that may have outlived its problem:

.. code-block:: python

   scalar_ff = phys.GetScalarFluxFieldFunction()[0]

   solver.Advance()
   if scalar_ff.CanUpdate():
       scalar_ff.Update()

Calling ``Update()`` refreshes the same field-function object. Interpolators
and export calls that already reference that object will see the refreshed data
the next time they execute.

This is mainly useful in timestep loops or repeated-solve workflows where
reusing the same interpolation or export setup is clearer than reconstructing
field functions each time.

Scalar Flux Field Functions
===========================

``GetScalarFluxFieldFunction``
------------------------------

Use :py:meth:`pyopensn.solver.LBSProblem.GetScalarFluxFieldFunction` to create
scalar-flux or flux-moment field functions from the current scalar-flux state.

The common case is scalar flux only:

.. code-block:: python

   scalar_ffs = phys.GetScalarFluxFieldFunction()

This returns one field function per energy group, each representing the zeroth
moment of the flux.

If higher moments are needed too:

.. code-block:: python

   ff_by_group_and_moment = phys.GetScalarFluxFieldFunction(
       only_scalar_flux=False,
   )

In that form:

* ``result[g][m]`` is the field function for group ``g`` and moment ``m``

.. note::

   Most visualization and response workflows begin with scalar flux. Higher
   moments are mainly for diagnostics and specialized analysis.

Derived Field Functions
=======================

``CreateFieldFunction``
-----------------------

Use :py:meth:`pyopensn.solver.LBSProblem.CreateFieldFunction` to create a named
scalar field function derived from a 1D cross section or from the special case
``"power"``.

Examples:

.. code-block:: python

   fission_rate_ff = phys.CreateFieldFunction("fission_rate", "sigma_f")
   power_ff = phys.CreateFieldFunction("power_generation", "power")

For a built-in or custom 1D XS name, this creates the field:

.. code-block:: text

   sum_g xs[g] * phi_g

at each spatial point.

If a power-normalized view is needed:

.. code-block:: python

   power_ff = phys.CreateFieldFunction(
       "power_generation_norm",
       "power",
       power_normalization_target=1.0,
   )

Important points:

* ``name`` is the field-function name assigned to the returned object
* ``xs_name`` can be a built-in 1D XS name, a custom XS name, or ``"power"``
* ``power_normalization_target`` is optional and affects only the returned
  field function

.. note::

   ``power_normalization_target`` applies a power-based scaling to the returned
   field function only. It does not rescale the solver's internal flux vectors,
   and it does not change field functions created earlier. If you want a
   normalized power field, use
   ``CreateFieldFunction("name", "power", power_normalization_target=...)``.
   The same argument can also be used for other derived fields such as
   ``sigma_f * phi`` when a power-normalized post-processed view is desired.

Angular Flux Field Functions
============================

``GetAngularFieldFunctionList``
-------------------------------

Use
:py:meth:`pyopensn.solver.DiscreteOrdinatesProblem.GetAngularFieldFunctionList`
to create field functions for selected angular-flux components.

Example:

.. code-block:: python

   ang_ffs = phys.GetAngularFieldFunctionList(groups=[0], angles=[0])

Important requirements:

* this is available on :py:class:`pyopensn.solver.DiscreteOrdinatesProblem`
* angular-flux storage must be enabled with ``save_angular_flux=True``
* the returned field functions are snapshots created from the currently stored
  angular flux
* if the stored angular flux changes later, call ``Update()`` on the existing
  angular field functions or create new ones from the current state

.. important::

   For transient problems, ``save_angular_flux=True`` is not optional. It is
   required by the transient solver itself, and it is also required if you
   want to create angular-flux field functions or query ``GetPsi()`` during the
   transient run.

This is mainly useful for:

* angular diagnostics,
* leakage studies,
* checking directional structure in difficult problems.

.. note::

   Angular flux is much larger than scalar flux. Only enable angular-flux
   storage and create angular field functions when the workflow actually needs
   them.

Scalar-flux field functions use names based on group and moment, for example:

* ``phi_g000_m00``
* ``phi_g001_m00``
* ``phi_g000_m01``

If the problem uses a field-function prefix, that prefix is applied to these
names as well.

Exporting Field Functions
=========================

The main Python export routine is:

* :py:meth:`pyopensn.fieldfunc.FieldFunctionGridBased.ExportMultipleToPVTU`

Export scalar flux:

.. code-block:: python

   from pyopensn.fieldfunc import FieldFunctionGridBased

   scalar_ffs = phys.GetScalarFluxFieldFunction()
   FieldFunctionGridBased.ExportMultipleToPVTU(
       scalar_ffs,
       "scalar_flux",
   )

Export a derived field:

.. code-block:: python

   from pyopensn.fieldfunc import FieldFunctionGridBased

   power_ff = phys.CreateFieldFunction("power_generation", "power")
   FieldFunctionGridBased.ExportMultipleToPVTU(
       [power_ff],
       "power",
   )

Export several fields together:

.. code-block:: python

   from pyopensn.fieldfunc import FieldFunctionGridBased

   scalar_ffs = phys.GetScalarFluxFieldFunction()
   power_ff = phys.CreateFieldFunction("power_generation", "power")

   FieldFunctionGridBased.ExportMultipleToPVTU(
       scalar_ffs + [power_ff],
       "transport_outputs",
   )

.. note::

   Export after the solve, not before. For transient problems, export after the
   completed timestep whose state you want to visualize.

Field-Function Interpolation
============================

Point Interpolation
-------------------

Use :py:class:`pyopensn.fieldfunc.FieldFunctionInterpolationPoint` to evaluate a
field function at a point.

Example:

.. code-block:: python

   from pyopensn.fieldfunc import FieldFunctionInterpolationPoint

   ffi = FieldFunctionInterpolationPoint()
   ffi.AddFieldFunction(phys.GetScalarFluxFieldFunction()[0])
   ffi.Initialize()
   ffi.Execute()
   value = ffi.GetPointValue()

This is useful for detector-like spot checks or comparison against an analytic
value.

Line Interpolation
------------------

Use :py:class:`pyopensn.fieldfunc.FieldFunctionInterpolationLine` to sample a
field function along a line.

Example:

.. code-block:: python

   from pyopensn.fieldfunc import FieldFunctionInterpolationLine
   from pyopensn.math import Vector3

   ffi = FieldFunctionInterpolationLine()
   ffi.AddFieldFunction(phys.GetScalarFluxFieldFunction()[0])
   ffi.SetInitialPoint(Vector3(0.0, 0.0, 0.0))
   ffi.SetFinalPoint(Vector3(10.0, 0.0, 0.0))
   ffi.SetNumberOfPoints(101)
   ffi.Initialize()
   ffi.Execute()
   ffi.ExportToCSV("centerline")

This is useful for profiles, attenuation curves, and one-dimensional
comparisons across runs.

Volume Interpolation
--------------------

Use :py:class:`pyopensn.fieldfunc.FieldFunctionInterpolationVolume` for region
averages, sums, maxima, and function-weighted variants over a logical volume.

Example:

.. code-block:: python

   from pyopensn.fieldfunc import FieldFunctionInterpolationVolume

   ffi = FieldFunctionInterpolationVolume()
   ffi.AddFieldFunction(phys.GetScalarFluxFieldFunction()[0])
   ffi.SetLogicalVolume(my_lv)
   ffi.SetOperationType("avg")
   ffi.Initialize()
   ffi.Execute()
   avg_value = ffi.GetValue()

Available operation types are:

* ``"sum"``
* ``"avg"``
* ``"max"``
* ``"sum_func"``
* ``"avg_func"``
* ``"max_func"``

The ``*_func`` variants use a scalar material function supplied with
``SetOperationFunction``.

Volume Postprocessor
====================

The :py:class:`pyopensn.post.VolumePostprocessor` computes scalar-flux integrals,
maxima, minima, or volume-weighted averages over spatial regions and energy
groups. It produces a single result value per region and group, making it useful
for reaction rates and monitoring quantities across geometry subsets or
energy ranges.

Basic Usage
-----------

Create and execute a volume postprocessor:

.. code-block:: python

   from pyopensn.post import VolumePostprocessor

   pps = VolumePostprocessor(
       problem=phys,
       value_type="integral",
   )
   pps.Execute()
   result = pps.GetValue()

Available operation types are:

* ``"integral"`` — volume-weighted integral of scalar flux
* ``"avg"`` — volume-weighted average of scalar flux
* ``"max"`` — maximum scalar flux in region
* ``"min"`` — minimum scalar flux in region

Spatial Restriction
-------------------

By default, a postprocessor operates over the entire domain. Restrict computation
to mesh blocks:

.. code-block:: python

   pps = VolumePostprocessor(
       problem=phys,
       value_type="integral",
       block_ids=[1, 2],
   )

Or restrict to one or more logical volumes:

.. code-block:: python

   from pyopensn.logvol import RPPLogicalVolume

   lv = RPPLogicalVolume(xmin=0.0, xmax=1.0, ymin=0.0, ymax=0.5)
   pps = VolumePostprocessor(
       problem=phys,
       value_type="avg",
       logical_volumes=[lv],
   )

Combine block and logical-volume restrictions — the postprocessor uses only cells
inside both:

.. code-block:: python

   lv1 = RPPLogicalVolume(xmin=0.0, xmax=1.0, ymin=0.0, ymax=0.5)
   lv2 = RPPLogicalVolume(xmin=1.0, xmax=2.0, ymin=0.0, ymax=0.5)
   pps = VolumePostprocessor(
       problem=phys,
       value_type="integral",
       logical_volumes=[lv1, lv2],
       block_ids=[1],
   )

Each logical volume produces one row of results.

Energy Restriction
-------------------

By default, a postprocessor returns results for all energy groups. Restrict to a
single group:

.. code-block:: python

   pps = VolumePostprocessor(
       problem=phys,
       value_type="integral",
       group=6,
   )

Or restrict to a single groupset:

.. code-block:: python

   pps = VolumePostprocessor(
       problem=phys,
       value_type="integral",
       groupset=1,
   )

The ``group`` and ``groupset`` parameters are mutually exclusive.

Multipliers and Cross-Section Weighting
----------------------------------------

By default, the postprocessor multiplies each group's scalar flux by 1.0. Apply a
uniform multiplier:

.. code-block:: python

   pps = VolumePostprocessor(
       problem=phys,
       value_type="integral",
       multiplier=2.5,
   )

Apply group-specific multipliers:

.. code-block:: python

   group_mults = [1.0, 1.5, 2.0, 2.0, 2.0, 2.0, 3.0, 3.0]  # one per group
   pps = VolumePostprocessor(
       problem=phys,
       value_type="integral",
       group_multipliers=group_mults,
   )

Weight by a cross section (for example, fission rate):

.. code-block:: python

   pps = VolumePostprocessor(
       problem=phys,
       value_type="integral",
       xs_multiplier="sigma_f",
   )

The cross-section name must exist in the problem's XS definitions. Only one of
``multiplier``, ``group_multipliers``, and ``xs_multiplier`` may be specified.

Results
-------

After calling ``Execute()``, retrieve results with ``GetValue()``. The return
value is a 2D array indexed as ``result[region][group]``:

.. code-block:: python

   pps = VolumePostprocessor(
       problem=phys,
       value_type="integral",
   )
   pps.Execute()
   result = pps.GetValue()

   # Single region, all groups
   # result[0] is a vector of values, one per group
   for group_index, value in enumerate(result[0]):
       print(f"Group {group_index}: {value}")

With multiple logical volumes:

.. code-block:: python

   lv1 = RPPLogicalVolume(...)
   lv2 = RPPLogicalVolume(...)
   pps = VolumePostprocessor(
       problem=phys,
       value_type="integral",
       logical_volumes=[lv1, lv2],
   )
   pps.Execute()
   result = pps.GetValue()

   # result[0] is values for lv1
   # result[1] is values for lv2

With energy restriction:

.. code-block:: python

   pps = VolumePostprocessor(
       problem=phys,
       value_type="integral",
       group=3,
   )
   pps.Execute()
   result = pps.GetValue()

   # result[0] is a vector with one element (the single group)
   print(result[0][0])

Other Useful Post-Processing Paths
==================================

For discrete ordinates problems, other useful output paths include:

* :py:meth:`GetPsi` for direct access to stored angular-flux arrays
* :py:meth:`ComputeLeakage` for boundary leakage
* :py:meth:`WriteAngularFluxes` and :py:meth:`ReadAngularFluxes` for angular
  flux storage and restart-style workflows

These are not field-function exports, but they are often part of the same
post-processing workflow.

Transient Workflows
===================

For transient problems, the usual pattern is:

1. advance the timestep,
2. update or create the desired field functions,
3. export or interpolate them,
4. repeat as needed.

If the transient workflow needs angular-flux output as well, the problem must
have been created with ``options={"save_angular_flux": True}`` from the start.

For example:

.. code-block:: python

   scalar_ff = phys.GetScalarFluxFieldFunction()[0]

   while not solver.Finished():
       solver.Advance()
       scalar_ff.Update()
       ...

This keeps one field-function object and refreshes it from each completed
timestep. Creating a new field function after each ``Advance()`` is also valid,
but is usually unnecessary when the existing object supports ``Update()``.

Practical Guidance
==================

For most workflows:

* use :py:meth:`GetScalarFluxFieldFunction` for scalar flux,
* use :py:meth:`CreateFieldFunction` for power or XS-weighted derived outputs,
* use :py:meth:`GetAngularFieldFunctionList` only when angular information is
  actually needed,
* use :py:meth:`ExportMultipleToPVTU` for visualization output,
* use point interpolation for spot checks,
* use line interpolation for profiles,
* use volume interpolation for averages and integrated responses.
* call ``Update()`` before reusing a field function after the solver state
  changes, or create a fresh field function from the current state.

If a script is becoming complicated, it is often worth separating the solve and
the post-processing logic into different helper functions. That keeps the
transport setup readable and makes the output workflow easier to reuse.
