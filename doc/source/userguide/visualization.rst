=============
Visualization
=============

Visualization in OpenSn usually means exporting field functions or meshes to
VTK-family files and then opening those files in a compatible viewer (VisIt,
ParaView, etc.).

This section is intentionally brief because most of the mechanics are covered by
the field-function and mesh APIs themselves. The aim here is to show a typical
visualization workflow.

Overview
========

The OpenSn export methods are:

* :py:meth:`pyopensn.fieldfunc.FieldFunctionGridBased.ExportMultipleToPVTU`
* :py:meth:`pyopensn.mesh.MeshContinuum.ExportToPVTU`

In practice, users usually visualize one of two things:

* the mesh itself, to verify geometry, block ids, and partitions,
* field functions, to inspect scalar flux, power, or other solution outputs.

Visualizing the Mesh
====================

Before solving a transport problem, it is often useful to export the mesh:

.. code-block:: python

   mesh.ExportToPVTU("mesh_check")

This is especially useful for:

* verifying imported geometry,
* checking block id assignment,
* checking partitioning on a parallel run.

.. note::

   Visualizing the mesh before solving is often the fastest way to catch setup
   mistakes. If the geometry or labels are wrong, the transport results will not
   become more trustworthy later.

Visualizing Scalar Flux
=======================

The most common field-function export is scalar flux:

.. code-block:: python

   from pyopensn.fieldfunc import FieldFunctionGridBased

   scalar_ffs = phys.GetScalarFluxFieldFunction()
   FieldFunctionGridBased.ExportMultipleToPVTU(
       scalar_ffs,
       "scalar_flux",
   )

This produces a parallel VTU dataset with each group's scalar flux written as a
field in the output set.

Visualizing Power
=================

Power can also be exported by creating a power field function on demand:

.. code-block:: python

   from pyopensn.fieldfunc import FieldFunctionGridBased

   power_ff = phys.CreateFieldFunction("power_generation", "power")
   FieldFunctionGridBased.ExportMultipleToPVTU(
       [power_ff],
       "power",
   )

This is commonly used in k-eigenvalue workflows or any problem where a
power-like field is part of the desired output.

Exporting Several Fields Together
=================================

It is often convenient to export several related field functions in one file
set:

.. code-block:: python

   from pyopensn.fieldfunc import FieldFunctionGridBased

   scalar_ffs = phys.GetScalarFluxFieldFunction()
   power_ff = phys.CreateFieldFunction("power_generation", "power")

   FieldFunctionGridBased.ExportMultipleToPVTU(
       scalar_ffs + [power_ff],
       "transport_outputs",
   )

This keeps the outputs together in one visualization dataset.

Transient Visualization
=======================

For transient problems, the usual pattern is:

1. advance the timestep,
2. update or retrieve the desired field functions,
3. export them,
4. repeat as needed.

This is easy to do in a manual Python loop, for example when only selected
time steps should be exported.

For example, create the field functions once and refresh them after each
timestep:

.. code-block:: python

   from pyopensn.fieldfunc import FieldFunctionGridBased

   scalar_ffs = phys.GetScalarFluxFieldFunction()

   for step in range(num_steps):
       solver.Advance()
       for ff in scalar_ffs:
           ff.Update()
       FieldFunctionGridBased.ExportMultipleToPVTU(
           scalar_ffs,
           f"scalar_flux_{step:04d}",
       )

Creating new field functions after each timestep is also valid, but is usually
not necessary when the existing objects support ``Update()``.

Practical Advice
================

* Export the mesh first when debugging geometry or labeling.
* Export scalar flux first when debugging transport behavior.
* Add power export only when power is actually a useful output for the problem.
* Keep visualization output separate from the core solver setup so the input is
  easier to read and maintain.

.. note::

   Visualization is usually the last step in the workflow. If the solve or the
   setup is not trustworthy yet, more exported files will not fix that. First
   verify the geometry, materials, sources, boundaries, and convergence. Then
   generate plots.
