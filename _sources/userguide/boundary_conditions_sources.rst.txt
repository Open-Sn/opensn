===============================
Boundary Conditions and Sources
===============================

Boundary conditions and sources are the primary drivers of a transport problem:

* boundary conditions drive particles through domain faces
* point sources inject particles at a point in space
* volumetric sources inject particles over cells or logical volumes

.. note::

   All of these source types are isotropic in angle. Angularly structured inflow
   is handled through boundary conditions.

This section covers:

* the supported boundary-condition types
* point and volumetric source objects
* how these inputs are specified at construction time
* how to replace them later with the problem-level setter methods
* how time dependence works for transient problems

Overview
========

OpenSn handles boundaries and sources through the problem constructor and the
corresponding problem-level setter methods. The following interfaces are
available for specifying boundary conditions and sources:

* :py:meth:`pyopensn.solver.DiscreteOrdinatesProblem.SetBoundaryOptions`
* :py:meth:`pyopensn.solver.LBSProblem.SetPointSources`
* :py:meth:`pyopensn.solver.LBSProblem.SetVolumetricSources`
* :py:class:`pyopensn.source.PointSource`
* :py:class:`pyopensn.source.VolumetricSource`

In most inputs, the workflow is simple:

1. pass ``boundary_conditions``, ``point_sources``, and
   ``volumetric_sources`` in the problem constructor, or
2. create the problem first and update those inputs later with
   :py:meth:`SetBoundaryOptions`, :py:meth:`SetPointSources`, and
   :py:meth:`SetVolumetricSources`

For most users, the practical order is also straightforward:

1. start with simple boundary conditions such as ``"vacuum"``,
   ``"reflecting"``, or ``"isotropic"``,
2. add point or volumetric sources if the problem is source driven,
3. only use :py:class:`pyopensn.math.AngularFluxFunction` or
   :py:class:`pyopensn.math.VectorSpatialFunction` when the boundary or source
   really needs callback-based custom behavior.

.. note::

   ``AngularFluxFunction`` and ``VectorSpatialFunction`` are the flexible
   callback-based tools for arbitrary boundary conditions, but they are not the
   starting point for most inputs. Most problems use the built-in
   boundary-condition types and simple point or volumetric source objects.

Boundary Conditions
===================

Boundary conditions are defined on
:py:class:`pyopensn.solver.DiscreteOrdinatesProblem` using the
``boundary_conditions`` constructor argument or
:py:meth:`pyopensn.solver.DiscreteOrdinatesProblem.SetBoundaryOptions`.

Supported boundary types are:

* ``"vacuum"`` allows particles to leave the domain but does not prescribe any
  incoming flux
* ``"reflecting"`` returns outgoing particles back into the domain with the
  appropriate reflected direction
* ``"isotropic"`` prescribes a group-wise isotropic incoming flux on the
  boundary
* ``"arbitrary"`` prescribes incoming angular flux through a user callback, so
  the inflow can depend on group and direction

Each boundary entry is a dictionary with:

* ``name``: the boundary id on the mesh
* ``type``: one of the supported boundary-condition types
* ``group_strength``: required for ``"isotropic"``
* ``function``: required for ``"arbitrary"``

.. note::

   OpenSn normalizes quadrature weights so that the full quadrature set sums to
   ``1.0``. For boundary inputs, that means ``group_strength`` and arbitrary
   boundary callback return values should be interpreted directly in OpenSn's
   discrete angular-flux convention. Users should not add extra ``4*pi``,
   ``2*pi``, or weight-based scaling to compensate for a different quadrature
   normalization.

Example:

.. code-block:: python

   boundary_conditions = [
       {"name": "xmin", "type": "vacuum"},
       {"name": "xmax", "type": "reflecting"},
       {"name": "ymin", "type": "isotropic", "group_strength": [1.0, 0.0]},
   ]

Vacuum Boundaries
-----------------

``"vacuum"`` means no incoming angular flux is prescribed on that boundary.

Example:

.. code-block:: python

   {"name": "zmin", "type": "vacuum"}

Use vacuum boundaries when particles leaving the domain should not return.

.. note::

   Vacuum is the default physical choice for many shielding and source-driven
   transport problems.

Reflecting Boundaries
---------------------

``"reflecting"`` mirrors the outgoing angular flux back into the domain.

Example:

.. code-block:: python

   {"name": "xmin", "type": "reflecting"}

Use reflecting boundaries for symmetry planes or to model an infinite medium.

.. note::

   Reflecting boundaries are a modeling statement, not just a numerical trick.
   They are appropriate when the physical problem really has mirror symmetry or
   an intended periodic-like repeated structure.

Isotropic Boundaries
--------------------

``"isotropic"`` applies a group-wise isotropic incoming angular flux on the
boundary.

Example:

.. code-block:: python

   bsrc = [0.0] * num_groups
   bsrc[0] = 1.0

   {"name": "outside", "type": "isotropic", "group_strength": bsrc}

Important points:

* ``group_strength`` must have one entry per energy group
* the values are incoming boundary angular flux values
* because OpenSn normalizes quadrature weights to sum to ``1.0``, a value
  ``group_strength[g] = X`` means "apply a uniform isotropic incoming angular
  flux of value ``X`` in group ``g``" in OpenSn's discrete convention

.. note::

   Isotropic boundaries are often the simplest way to drive a benchmark or test
   problem when the intent is "inject a known incoming flux on this face" rather
   than "create particles throughout a volume."

Arbitrary Boundaries
--------------------

``"arbitrary"`` uses a Python callback wrapped in
:py:class:`pyopensn.math.AngularFluxFunction`.

The callback signature is:

.. code-block:: python

   def bc_func(group_index, direction_index) -> float:
       ...

The first argument is the global energy-group index. The second is the angular
quadrature direction index within the groupset quadrature.

The callback should return the incoming angular-flux value for that
group/direction pair. Because OpenSn normalizes the quadrature weights to sum
to ``1.0``, the return value should be given directly in OpenSn's discrete
angular-flux convention rather than rescaled by external solid-angle factors.

Example:

.. code-block:: python

   from pyopensn.math import AngularFluxFunction

   def xmin_bc_func(group_index, direction_index):
       if group_index == 0 and direction_index < 4:
           return 0.5
       return 0.0

   xmin_bc = AngularFluxFunction(xmin_bc_func)

   boundary_conditions = [
       {"name": "xmin", "type": "arbitrary", "function": xmin_bc},
       {"name": "xmax", "type": "vacuum"},
   ]

.. note::

   ``AngularFluxFunction`` is the most flexible boundary option, but it is also
   the easiest one to misuse. A good first check is always whether the callback
   is being written in terms of the correct quadrature direction indices.

.. note::

   Use ``"arbitrary"`` when the inflow really depends on direction. If the
   intended condition is just a group-wise isotropic inflow, the
   ``"isotropic"`` path is clearer and easier to review.

.. note::

   If you want an arbitrary boundary to behave like a uniform isotropic inflow
   of value ``X`` in one group, return ``X`` for every incoming direction in
   that group. Do not multiply by quadrature weights in the callback.

Setting and Replacing Boundary Conditions
-----------------------------------------

Boundary conditions can be provided in the problem constructor:

.. code-block:: python

   phys = DiscreteOrdinatesProblem(
       ...,
       boundary_conditions=[
           {"name": "xmin", "type": "vacuum"},
           {"name": "xmax", "type": "reflecting"},
       ],
   )

or updated later with :py:meth:`SetBoundaryOptions`:

.. code-block:: python

   phys.SetBoundaryOptions(
       clear_boundary_conditions=True,
       boundary_conditions=[
           {"name": "xmin", "type": "vacuum"},
           {"name": "xmax", "type": "vacuum"},
       ],
   )

Important behavior:

* ``clear_boundary_conditions=True`` removes the current set before applying the
  new list
* :py:meth:`LBSProblem.SetAdjoint` clears all boundary conditions during a mode
  transition, so they must be reapplied afterward

.. note::

   Boundary-condition updates are a natural way to steer a problem between
   solves. They are often clearer than rebuilding the entire problem object when
   only the boundary driving conditions are changing.

Point Sources
=============

:py:class:`pyopensn.source.PointSource` represents a multi-group isotropic point
source.

Constructor options:

* ``location``: required XYZ point location
* ``strength``: optional group-wise strength vector
* ``strength_function``: optional callable
* ``start_time``: optional activation time
* ``end_time``: optional deactivation time

Exactly one of ``strength`` or ``strength_function`` must be provided.

Constant Point Source
---------------------

Example:

.. code-block:: python

   from pyopensn.source import PointSource

   psrc = PointSource(
       location=(0.0, 0.0, 0.0),
       strength=[1.0, 0.0, 0.0],
   )

This defines an isotropic point source at the given coordinates with one value
per energy group.

.. note::

   Point sources are best used when the physical source is truly localized. If
   the source occupies a region of space, a volumetric source is usually the
   clearer model.

Time-Windowed Point Sources
---------------------------

For non-callback point sources, activity can be limited in time with
``start_time`` and ``end_time``:

.. code-block:: python

   psrc = PointSource(
       location=(0.0, 0.0, 0.0),
       strength=[1.0],
       start_time=0.0,
       end_time=1.0,
   )

The source is active when:

* ``time >= start_time``
* ``time <= end_time``

.. note::

   The start/end window is the simplest transient source option when the source
   turns on and off at known times but does not otherwise change shape.

Functional Point Sources
------------------------

Instead of a fixed vector, a point source can use ``strength_function``.

The callback may accept either:

* ``(group, time)`` for transient problems, or
* ``(group)`` for steady-state use

Example:

.. code-block:: python

   def point_strength(group, time):
       if group == 0:
           return 1.0 + time
       return 0.0

   psrc = PointSource(
       location=(0.0, 0.0, 0.0),
       strength_function=point_strength,
   )

Important restriction:

* if ``strength_function`` is provided, do not also provide ``start_time`` or
  ``end_time``; the time dependence should be handled inside the callback

.. note::

   A callback-based point source is the right choice when the source strength is
   genuinely time-dependent or needs Python-side logic. For simple on/off
   behavior, the fixed-strength plus start/end-time path is easier to read.

Volumetric Sources
==================

:py:class:`pyopensn.source.VolumetricSource` represents a multi-group isotropic
volumetric source.

It can be defined over:

* one or more ``block_ids``
* a ``logical_volume``
* or the intersection of both

Constructor options:

* ``block_ids``
* ``logical_volume``
* ``group_strength``
* ``func``
* ``strength_function``
* ``start_time``
* ``end_time``

A volumetric source must specify at least one spatial selector:

* ``block_ids``, or
* ``logical_volume``

and exactly one source-definition mode:

* ``group_strength``
* ``func``
* ``strength_function``

Block-ID-Based Volumetric Sources
---------------------------------

Example:

.. code-block:: python

   from pyopensn.source import VolumetricSource

   src = VolumetricSource(
       block_ids=[0, 1],
       group_strength=[5.0, 0.0],
   )

This applies a constant isotropic volumetric source in all local cells whose
block id is 0 or 1.

.. note::

   Block-id-based volumetric sources are usually the cleanest option when the
   source region already matches the material or mesh-region labeling.

Logical-Volume-Based Volumetric Sources
---------------------------------------

Example:

.. code-block:: python

   from pyopensn.logvol import RPPLogicalVolume
   from pyopensn.source import VolumetricSource

   source_region = RPPLogicalVolume(
       xmin=0.0, xmax=1.0,
       ymin=0.0, ymax=1.0,
       zmin=0.0, zmax=1.0,
   )

   src = VolumetricSource(
       logical_volume=source_region,
       group_strength=[1.0],
   )

This is useful when the source region is geometric and should not depend only on
the block id layout.

If both ``logical_volume`` and ``block_ids`` are supplied, the source acts only
in cells that satisfy both conditions.

.. note::

   The combined ``logical_volume`` plus ``block_ids`` path is useful for
   selecting a subset of a material region without changing the underlying mesh
   labeling.

Spatial Functional Volumetric Sources
-------------------------------------

Use :py:class:`pyopensn.math.VectorSpatialFunction` when the source varies in
space.

The callback signature is:

.. code-block:: python

   def spatial_source(point, num_groups) -> list[float]:
       ...

Example:

.. code-block:: python

   from pyopensn.math import VectorSpatialFunction

   def source_profile(point, num_groups):
       value = point.x + point.y
       return [value] * num_groups

   src = VolumetricSource(
       block_ids=[0],
       func=VectorSpatialFunction(source_profile),
   )

Important behavior:

* the function returns one value per group
* the spatial function is evaluated only in the selected cells
* a spatial function is mutually exclusive with ``group_strength`` and
  ``strength_function``

.. note::

   ``VectorSpatialFunction`` is the natural choice when the source has a spatial
   profile but does not need custom angular structure.

Group-Time Functional Volumetric Sources
----------------------------------------

Use ``strength_function`` when the volumetric source varies by group and time
but not by space within its selected region.

The callback may accept either:

* ``(group, time)`` for transient problems, or
* ``(group)`` for steady-state use

Example:

.. code-block:: python

   def group_time_strength(group, time):
       if group == 0:
           return 1.0 if time < 0.5 else 0.0
       return 0.0

   src = VolumetricSource(
       block_ids=[0],
       strength_function=group_time_strength,
   )

As with point sources:

* do not combine ``strength_function`` with ``start_time`` or ``end_time``

.. note::

   ``strength_function`` is best when the source region is fixed but its
   magnitude changes over time. If both the region and the profile vary, update
   the source object itself between solves or timesteps instead.

Time-Windowed Volumetric Sources
--------------------------------

For constant-strength or spatial-function sources, ``start_time`` and
``end_time`` provide a simple activation window.

Example:

.. code-block:: python

   src = VolumetricSource(
       block_ids=[0],
       group_strength=[1.5],
       start_time=0.0,
       end_time=1.0,
   )

This is a common pattern in transient tests.

.. note::

   When the source should simply be present during a known time interval,
   ``start_time`` and ``end_time`` are easier to understand and review than a
   callback that reproduces the same logic.

Attaching Sources at Problem Construction
=========================================

Both point and volumetric sources can be passed directly to the problem
constructor.

Example:

.. code-block:: python

   psrc = PointSource(location=(0.0, 0.0, 0.0), strength=[1.0])
   vsrc = VolumetricSource(block_ids=[0], group_strength=[2.0])

   phys = DiscreteOrdinatesProblem(
       ...,
       point_sources=[psrc],
       volumetric_sources=[vsrc],
       boundary_conditions=[
           {"name": "zmin", "type": "vacuum"},
           {"name": "zmax", "type": "vacuum"},
       ],
   )

This is the cleanest pattern when the driving terms are known up front.

Replacing Sources After Construction
====================================

Problem objects also expose:

* :py:meth:`pyopensn.solver.LBSProblem.SetPointSources`
* :py:meth:`pyopensn.solver.LBSProblem.SetVolumetricSources`

These methods accept:

* ``clear_point_sources=True`` or ``clear_volumetric_sources=True``
* a new source list

Examples:

.. code-block:: python

   phys.SetPointSources(
       clear_point_sources=True,
       point_sources=[psrc_new],
   )

   phys.SetVolumetricSources(
       clear_volumetric_sources=True,
       volumetric_sources=[vsrc_new],
   )

.. note::

   Updating sources through the problem object is often the simplest way to
   drive a sequence of solves or a Python-controlled transient without
   reconstructing the entire problem.

Transient Source Behavior
=========================

For transient problems, time dependence can enter through:

* ``start_time`` / ``end_time``
* ``strength_function`` callbacks
* explicit replacement of source objects between calls to
  :py:meth:`Advance`

Example of a time-windowed volumetric source:

.. code-block:: python

   src = VolumetricSource(
       block_ids=[0],
       group_strength=[1.2],
       start_time=0.0,
       end_time=10.0,
   )

Example of replacing sources in a Python timestep loop:

.. code-block:: python

   solver.Initialize()

   for step in range(num_steps):
       if step == 10:
           phys.SetVolumetricSources(
               clear_volumetric_sources=True,
               volumetric_sources=[new_source],
           )
       solver.Advance()

.. note::

   There are two different design styles for transient sources. One is to keep a
   single source object and let it handle its own time dependence. The other is
   to replace source objects explicitly in the Python loop. Both are valid; the
   clearer choice depends on how much logic the source behavior really needs.

Boundary Names and Mesh Labels
==============================

Boundary-condition entries refer to mesh boundary ids, not geometric
directions in the abstract.

Typical boundary ids on orthogonal Cartesian meshes are:

* ``xmin``, ``xmax``
* ``ymin``, ``ymax``
* ``zmin``, ``zmax``

But meshes can also use custom boundary ids, for example:

* ``outside``
* ``left``
* ``right``

The boundary ids used in the input must match the boundary ids present on the
mesh.

Best Practices
==============

* Use vacuum and reflecting boundaries whenever they represent the actual
  physical model; they are the clearest and least error-prone options.
* Use isotropic boundaries for simple incoming-flux problems.
* Use arbitrary boundaries only when the inflow genuinely depends on angle.
* Use volumetric sources based on block ids when source regions already align
  with the mesh labeling.
* Use logical volumes when the source region is geometric rather than material
  based.
* Prefer ``start_time``/``end_time`` for simple on/off source behavior.
* Prefer callbacks when source strength must vary continuously with group or
  time.
* After calling :py:meth:`LBSProblem.SetAdjoint`, reapply sources
  and boundary conditions before solving again.

.. note::

   The cleanest source and boundary inputs are usually the ones that say exactly
   what the physics is with the least machinery. Start simple, then move to
   callbacks or explicit update logic only when the problem actually needs that
   flexibility.
