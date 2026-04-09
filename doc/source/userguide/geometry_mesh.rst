=================
Geometry and Mesh
=================

In OpenSn, geometry is represented through the mesh. There is no separate
CAD-style geometry object that the transport solver uses directly. Instead,
users create or import a mesh with block ids and boundary ids, and then refer
to those ids when defining materials, sources, and boundary conditions.

Logical volumes provide an additional geometry-like tool. They are especially
useful when the mesh already exists and regions still need to be labeled for
materials, boundaries, sources, or post-processors.

This section explains the main mesh objects, supported mesh generators,
partitioning choices, and the workflows that matter most to transport problems.

Overview
========

The Python geometry interface API is:

* :py:class:`pyopensn.mesh.MeshContinuum`
* :py:class:`pyopensn.mesh.MeshGenerator`
* :py:class:`pyopensn.mesh.OrthogonalMeshGenerator`
* :py:class:`pyopensn.mesh.FromFileMeshGenerator`
* :py:class:`pyopensn.mesh.ExtruderMeshGenerator`
* :py:class:`pyopensn.mesh.SplitFileMeshGenerator`
* :py:class:`pyopensn.mesh.DistributedMeshGenerator`
* :py:class:`pyopensn.mesh.KBAGraphPartitioner`
* :py:class:`pyopensn.mesh.LinearGraphPartitioner`
* :py:class:`pyopensn.mesh.PETScGraphPartitioner`
* :py:class:`pyopensn.logvol.RPPLogicalVolume`
* :py:class:`pyopensn.logvol.RCCLogicalVolume`
* :py:class:`pyopensn.logvol.SphereLogicalVolume`
* :py:class:`pyopensn.logvol.BooleanLogicalVolume`
* :py:class:`pyopensn.logvol.SurfaceMeshLogicalVolume`

The most common workflow is:

1. create a mesh generator,
2. execute it to obtain a :py:class:`MeshContinuum`,
3. assign block ids and boundary ids if not already present on the mesh,
4. attach materials and sources using those ids,
5. build the transport problem on that mesh.

.. note::

   The mesh is not just a geometric container. It is also where the labeling
   information lives that later drives ``xs_map``, boundary conditions, and many
   source definitions. That is why mesh setup deserves its own careful pass in
   a new input.

MeshContinuum
=============

:py:class:`pyopensn.mesh.MeshContinuum` is the mesh object used by the transport
solvers. ``MeshContinuum`` is produced by a mesh generator and consumed by the
transport problem. It has a number of methods available for manipulating block
and boundary ids, querying properties, and exporting to file:

* :py:meth:`GetGlobalNumberOfCells`
* :py:meth:`SetUniformBlockID`
* :py:meth:`SetBlockIDFromLogicalVolume`
* :py:meth:`SetBlockIDFromFunction`
* :py:meth:`SetUniformBoundaryID`
* :py:meth:`SetBoundaryIDFromLogicalVolume`
* :py:meth:`SetOrthogonalBoundaries`
* :py:meth:`ComputeVolumePerBlockID`
* :py:meth:`ExportToPVTU`

.. note::

   A surprising amount of later trouble comes from mesh labeling rather than the
   transport algorithm. It is worth verifying block ids and boundary ids early
   before moving on to materials and groupsets.

Mesh Generators
===============

All mesh generators derive from :py:class:`pyopensn.mesh.MeshGenerator` and are
executed with:

.. code-block:: python

   mesh = generator.Execute()

Several mesh generators accept an ``inputs`` list of other mesh generators.
This allows generator pipelines where one generator consumes the output of
another.

This is the key design pattern behind mesh construction in the Python API:

* mesh generators are composable,
* generator execution is explicit,
* the transport problem only sees the final ``MeshContinuum``.

For example, these are all valid chaining patterns:

.. code-block:: python

   # Import an existing mesh directly.
   mesh = FromFileMeshGenerator(filename="mesh.msh").Execute()

   # Import a 2D mesh and then extrude it to 3D.
   base = FromFileMeshGenerator(filename="base_2d.obj")
   mesh = ExtruderMeshGenerator(
       inputs=[base],
       layers=[{"z": 1.0, "n": 4}],
   ).Execute()

   # Build an orthogonal mesh and then distribute it across ranks.
   mesh = DistributedMeshGenerator(
       inputs=[
           OrthogonalMeshGenerator(node_sets=[x_nodes, y_nodes, z_nodes]),
       ],
   ).Execute()

   # Build an orthogonal mesh and then write or read it through the split-file path.
   mesh = SplitFileMeshGenerator(
       inputs=[
           OrthogonalMeshGenerator(node_sets=[x_nodes, y_nodes, z_nodes]),
       ],
   ).Execute()

.. note::

   Chaining is one of the most useful parts of the mesh API. It lets you keep
   each step focused: one generator can create or import the base mesh, another
   can transform it, and another can handle partitioning or distribution.

The generators below are the main paths by which a ``MeshContinuum`` is created.

OrthogonalMeshGenerator
-----------------------

Use :py:class:`pyopensn.mesh.OrthogonalMeshGenerator` to build a structured
orthogonal mesh directly in Python.

Important constructor parameters:

* ``node_sets``: required list of coordinate arrays
* ``coord_sys``: ``"cartesian"`` or ``"cylindrical"``
* ``inputs``: optional list of preceding mesh generators
* ``partitioner``: optional graph partitioner
* ``replicated_mesh``: optional replicated/distributed choice

Example:

.. code-block:: python

   from pyopensn.mesh import OrthogonalMeshGenerator

   mesh = OrthogonalMeshGenerator(
       node_sets=[
           [0.0, 0.5, 1.0, 2.0],
           [0.0, 1.0, 2.0],
           [0.0, 1.0],
       ],
   ).Execute()

For orthogonal meshes, it is necessary to assign the boundary and material ids
after creation:

.. code-block:: python

   mesh.SetOrthogonalBoundaries()
   mesh.SetUniformBlockID(N)

This creates the usual named boundaries such as ``xmin``, ``xmax``, ``ymin``,
``ymax``, ``zmin``, ``zmax``, and sets the block id of all cells to the integer
value N.

.. note::

   Orthogonal meshes are often the best way to start a new transport model.
   They keep the geometry simple, make the boundary naming obvious, and reduce
   the number of moving parts during early debugging.

FromFileMeshGenerator
---------------------

Use :py:class:`pyopensn.mesh.FromFileMeshGenerator` to import a mesh from an
external file or supported case directory. The mesh generator will automatically
determine the file type.

Important constructor parameters:

* ``filename``: required input path
* ``inputs``: optional list of preceding mesh generators
* ``partitioner``: optional graph partitioner
* ``replicated_mesh``: optional replicated/distributed choice

Supported import families include:

* ``.obj``: Wavefront OBJ exporters from Blender, MeshLab, FreeCAD, and many
  CAD or geometry tools.
* ``.msh``: Gmsh.
* ``.e``: ExodusII-producing tools such as Cubit, Coreform Cubit, SEACAS
  tools, and some finite-element preprocessors.
* ``.vtu``: VTK or ParaView writers, ``meshio``, and many FEM or CFD tools
  that export VTK Unstructured Grid files.
* ``.pvtu``: ParaView or VTK parallel XML outputs, usually produced by
  parallel VTK writers rather than authored directly.
* ``.case``: EnSight Gold case format exporters from ParaView, EnSight, and
  some simulation codes.
* OpenFOAM case directories with the expected ``constant/polyMesh`` content:
  OpenFOAM utilities such as ``blockMesh``, ``snappyHexMesh``, and conversion
  tools such as ``gmshToFoam``.

Example:

.. code-block:: python

   from pyopensn.mesh import FromFileMeshGenerator

   mesh = FromFileMeshGenerator(
       filename="mesh.msh",
   ).Execute()

.. note::

   When importing a mesh, do not assume the block ids and boundary ids already
   match the intended transport input. Imported geometry frequently needs a
   cleanup or relabeling step.

ExtruderMeshGenerator
---------------------

Use :py:class:`pyopensn.mesh.ExtruderMeshGenerator` to extrude a lower-
dimensional input mesh into a higher-dimensional mesh.

Important constructor parameters:

* ``layers``: required list of layer dictionaries
* ``inputs``: mesh-generator inputs
* ``partitioner``: optional graph partitioner
* ``replicated_mesh``: optional replicated/distributed choice

Each layer dictionary uses:

* ``n``: number of sublayers
* exactly one of ``h`` or ``z``

Here:

* ``h`` means the total thickness of that layer segment measured from the
  current extrusion front.
* ``z`` means the absolute top coordinate to extrude to for that layer
  segment.

So if the current front is at ``z = 1.0``:

* ``{"n": 2, "h": 0.5}`` creates a segment of total thickness ``0.5`` and
  ends at ``z = 1.5``.
* ``{"n": 2, "z": 1.5}`` creates a segment that ends at the absolute location
  ``z = 1.5``.

In both cases, ``n`` controls how many sublayers that segment is divided into.

Example:

.. code-block:: python

   from pyopensn.mesh import FromFileMeshGenerator, ExtruderMeshGenerator

   base = FromFileMeshGenerator(filename="base_2d.msh")
   mesh = ExtruderMeshGenerator(
       inputs=[base],
       layers=[
           {"n": 2, "h": 0.5},
           {"n": 4, "h": 1.0},
       ],
   ).Execute()

This is a practical choice when:

* the natural geometry description is 2D,
* but the transport problem needs a three-dimensional mesh,
* and the geometry is uniform in the extrusion direction.

In other words, the 2D base mesh must be cleanly extrudable into the third
dimension without geometric variation from layer to layer. Block ids or other
material ids may vary by layer, but the cross-sectional geometry itself
must remain the same along the extrusion direction.

SplitFileMeshGenerator
----------------------

Use :py:class:`pyopensn.mesh.SplitFileMeshGenerator` when the mesh has already
been partitioned and written to disk as one mesh file per partition or rank.
Each rank reads the mesh data for its partition and constructs only its local
portion in memory. This is useful for large parallel problems where the
mesh is already available in per-partition form, since it avoids the need for
rank 0 to read the full mesh and distribute it. This generator can also be used
to write per-partition mesh files from an input mesh.

Important constructor parameters:

* ``file_base``
* ``inputs``
* ``partitioner``
* ``replicated_mesh``

Example:

.. code-block:: python

   mesh = SplitFileMeshGenerator(
       inputs=[
           OrthogonalMeshGenerator(node_sets=[x_nodes, y_nodes, z_nodes]),
       ],
   ).Execute()

DistributedMeshGenerator
------------------------

Use :py:class:`pyopensn.mesh.DistributedMeshGenerator` when a large mesh should
be partitioned and distributed across MPI ranks during mesh generation rather
than built in full on every rank. Rank 0 reads the full mesh from the
preceding generator in the chain, partitions it, and distributes each partition
to the MPI rank that owns it. This reduces memory usage and avoids requiring
every rank to first read or build the full mesh before the final distributed
mesh is assembled.

Conceptually, this is the in-memory analog of
:py:class:`pyopensn.mesh.SplitFileMeshGenerator`. The split-file path starts
from mesh data that is already written out per partition, while the distributed
generator starts from a single in-memory mesh and performs the partitioning and
distribution internally.

Important constructor parameters:

* ``inputs``
* ``partitioner``
* ``replicated_mesh``

Example:

.. code-block:: python

   mesh = DistributedMeshGenerator(
       inputs=[
           OrthogonalMeshGenerator(node_sets=[x_nodes, y_nodes, z_nodes]),
       ],
   ).Execute()

This generator is mainly a scalability tool rather than the first mesh-
construction choice for new users.

Coordinate Systems
==================

Several mesh generators accept ``coord_sys``:

* ``"cartesian"``
* ``"cylindrical"``

Cartesian is the default for ordinary ``x``, ``y``, ``z`` geometry. Use
``"cylindrical"`` for axisymmetric ``r-z`` style meshes and problems.

Examples:

.. code-block:: python

   mesh = OrthogonalMeshGenerator(
       node_sets=[x_nodes, y_nodes, z_nodes],
       coord_sys="cartesian",
   ).Execute()

   mesh = OrthogonalMeshGenerator(
       node_sets=[r_nodes, z_nodes],
       coord_sys="cylindrical",
   ).Execute()

   mesh = FromFileMeshGenerator(
       filename="../../../../assets/mesh/rz_rect_single.msh",
       coord_sys="cylindrical",
   ).Execute()

Partitioning
============

OpenSn partitions meshes for parallel execution using a graph partitioner.

Available partitioners in the Python API are:

* :py:class:`pyopensn.mesh.PETScGraphPartitioner`
* :py:class:`pyopensn.mesh.KBAGraphPartitioner`
* :py:class:`pyopensn.mesh.LinearGraphPartitioner`

If no partitioner is supplied, the mesh generators default to a PETSc-based
``parmetis`` partitioner.

PETScGraphPartitioner
---------------------

This is the general-purpose graph partitioner invoked through PETSc.

Use this first for most imported or unstructured meshes.

The available Python-facing ``type`` values are:

* ``"parmetis"``
* ``"ptscotch"``

ParMETIS works from the cell-adjacency graph and tries to minimize edge cuts
while balancing work across partitions. That makes it a strong default for
general unstructured or imported meshes.

PT-Scotch serves a similar role through the Scotch family of graph
partitioners. In practice, it is another general-purpose graph-based choice
for problems where you want an alternative to ParMETIS.

Example:

.. code-block:: python

   mesh = FromFileMeshGenerator(
       filename="mesh.msh",
       partitioner=PETScGraphPartitioner(type="parmetis"),
   ).Execute()

   mesh = FromFileMeshGenerator(
       filename="mesh.msh",
       partitioner=PETScGraphPartitioner(type="ptscotch"),
   ).Execute()

KBAGraphPartitioner
-------------------

KBA (Koch-Baker-Alcouffe) is an overlayed orthogonal-grid-based partitioner.
This partitioner overlays a rectangular partitioning structure on the domain
using user-specified cell counts and cut locations in ``x``, ``y``, and ``z``.
It is most natural for orthogonal meshes and for meshes that come from
structured or extruded workflows where a logically rectangular decomposition
still makes sense.

This is often a very good choice for orthogonal meshes because the user can
control the partition layout directly.

Example:

.. code-block:: python

   mesh = OrthogonalMeshGenerator(
       node_sets=[x_nodes, y_nodes, z_nodes],
       partitioner=KBAGraphPartitioner(
           nx=2,
           ny=2,
           nz=2,
           xcuts=[0.0],
           ycuts=[0.0],
           zcuts=[1.1],
       ),
   ).Execute()

LinearGraphPartitioner
----------------------

This is mainly useful for simple workflows and testing.

It should not be the general default for large realistic problems.

It partitions cells by their linear global ordering rather than by analyzing
the mesh graph. That can be acceptable for some simple orthogonal meshes, but
it is usually a poor choice for unstructured problems.

Example:

.. code-block:: python

   mesh = FromFileMeshGenerator(
       filename="mesh.msh",
       partitioner=LinearGraphPartitioner(),
   ).Execute()

   # Or force everything to a single rank for testing.
   mesh = FromFileMeshGenerator(
       filename="mesh.msh",
       partitioner=LinearGraphPartitioner(all_to_rank=0),
   ).Execute()

.. note::

   Most users should not start with partitioner tuning. First get the problem
   running correctly on a single rank or with a simple partitioning choice. Then
   tune decomposition if the problem size justifies it.

Replicated Meshes
=================

The main mesh generators support ``replicated_mesh``.

When ``replicated_mesh=False``:

* the final :py:class:`MeshContinuum` on each rank contains only its locally
  owned cells plus the ghost cells it needs.

When ``replicated_mesh=True``:

* the full mesh is present on every rank.

This flag affects the final ``MeshContinuum``, not the earlier mesh-generation
workflow.

For the ordinary generator path, ranks still participate in creating or reading
the full intermediate mesh before the final ``MeshContinuum`` is built. The
difference is what survives into the final mesh object:

* with ``replicated_mesh=False``, each rank keeps only its local piece plus
  ghosts,
* with ``replicated_mesh=True``, each rank keeps the full mesh.

If you want to avoid having every rank read or build the full mesh in the first
place, use :py:class:`pyopensn.mesh.DistributedMeshGenerator` or
:py:class:`pyopensn.mesh.SplitFileMeshGenerator`.

Distributed meshes are the normal choice for large production problems.
Replicated meshes are often more convenient for:

* small problems,
* debugging,
* early development and geometry verification.

Block IDs
=========

Block ids are integer labels stored on cells. They are used most prominently by
the ``xs_map`` material-assignment system in transport problems.

For meshes imported with :py:class:`pyopensn.mesh.FromFileMeshGenerator`, block
ids should have been added by the external mesh-generation tool. The methods
below are primarily for use with orthogonal meshes constructed via the input.

The main block-id assignment methods are:

* :py:meth:`SetUniformBlockID`
* :py:meth:`SetBlockIDFromLogicalVolume`
* :py:meth:`SetBlockIDFromFunction`

Uniform assignment
------------------

Use ``SetUniformBlockID`` when the entire mesh belongs to one material region:

.. code-block:: python

   mesh.SetUniformBlockID(0)

Logical-volume assignment
-------------------------

Use ``SetBlockIDFromLogicalVolume`` when a geometric region should receive a
specific block id:

.. code-block:: python

   mesh.SetBlockIDFromLogicalVolume(fuel_lv, 0, True)

Function-based assignment
-------------------------

Use ``SetBlockIDFromFunction`` when the labeling rule is easiest to express as
a Python function of position.

.. note::

   For most user inputs, the important thing is not how the block ids are
   assigned but that they are assigned consistently and then used consistently
   in ``xs_map``.

Boundary Naming and Boundary IDs
================================

OpenSn transport inputs refer to boundaries by id, not raw geometric faces, so
meaningful boundary assignment matters.

For meshes imported with :py:class:`pyopensn.mesh.FromFileMeshGenerator`,
boundary ids should have been set by the external meshing tool. The methods
below are primarily for use with orthogonal meshes constructed via the input.

The main boundary-assignment methods are:

* :py:meth:`SetUniformBoundaryID`
* :py:meth:`SetBoundaryIDFromLogicalVolume`
* :py:meth:`SetOrthogonalBoundaries`

For orthogonal meshes, ``SetOrthogonalBoundaries`` is usually the simplest
choice.

For more general workflows, logical-volume boundary assignment is often the
most practical way to label selected surfaces.

.. note::

   Boundary id mistakes are common because they do not always look like
   geometry mistakes at first. A vacuum boundary on the wrong face can look like
   a physics problem when it is really just a naming mismatch.

Logical Volumes
===============

Logical volumes are geometric selectors. They answer the question:

"Is this point inside the region?"

A cell is considered inside a logical volume if the cell's centroid is inside
the volume.

Available logical-volume classes are:

* :py:class:`RPPLogicalVolume`
* :py:class:`RCCLogicalVolume`
* :py:class:`SphereLogicalVolume`
* :py:class:`BooleanLogicalVolume`
* :py:class:`SurfaceMeshLogicalVolume`

Logical volumes are used for:

* block-id assignment,
* boundary-id assignment,
* source-region selection,
* post-processing region selection.

This makes them one of the most important geometry tools outside the mesh
itself.

For mesh labeling, the key methods are:

* :py:meth:`MeshContinuum.SetBlockIDFromLogicalVolume`
* :py:meth:`MeshContinuum.SetBoundaryIDFromLogicalVolume`

In both cases, the ``inside`` argument controls whether the selected cells or
faces are the ones whose centroids are inside the logical volume or outside it.
This is useful when you want to define either the interior region directly or
its complement.

.. note::

   Logical-volume selection is centroid-based. That is usually what you want,
   but it is important to remember on coarse meshes or for very thin regions:
   a cell is not selected because it intersects the volume, only because its
   centroid lies inside it.

Basic usage pattern
-------------------

For block-id assignment:

.. code-block:: python

   from pyopensn.logvol import RPPLogicalVolume

   fuel_region = RPPLogicalVolume(
       xmin=-0.5, xmax=0.5,
       ymin=-0.5, ymax=0.5,
       infz=True,
   )
   mesh.SetBlockIDFromLogicalVolume(fuel_region, 1, True)

For boundary labeling:

.. code-block:: python

   left_face_region = RPPLogicalVolume(
       xmin=-1.0, xmax=-1.0,
       infy=True,
       infz=True,
   )
   mesh.SetBoundaryIDFromLogicalVolume(left_face_region, "xmin", True)

For selecting the complement of a region:

.. code-block:: python

   mesh.SetBlockIDFromLogicalVolume(fuel_region, 0, False)

Here, ``False`` means "apply this assignment to cells whose centroids are
outside the logical volume."

``RPPLogicalVolume``
--------------------

Use :py:class:`RPPLogicalVolume` for rectangular-parallelepiped style regions.
This is often the most convenient logical volume for axis-aligned box-like
regions, slabs, strips, and half-spaces.

Important constructor parameters:

* ``xmin``, ``xmax``
* ``ymin``, ``ymax``
* ``zmin``, ``zmax``
* ``infx``, ``infy``, ``infz``

The ``inf*`` flags let you ignore bounds in one or more directions. That makes
``RPPLogicalVolume`` useful not only for finite boxes but also for infinite
slabs and strips.

Examples:

.. code-block:: python

   # A finite box.
   box = RPPLogicalVolume(
       xmin=-0.5, xmax=0.5,
       ymin=-0.5, ymax=0.5,
       zmin=0.0, zmax=1.0,
   )

   # An infinite slab in y and z, bounded only in x.
   slab = RPPLogicalVolume(
       xmin=-0.25, xmax=0.25,
       infy=True,
       infz=True,
   )

   # A whole-domain selector.
   whole_domain = RPPLogicalVolume(infx=True, infy=True, infz=True)

Use this volume when the region is naturally described by coordinate-aligned
limits.

``RCCLogicalVolume``
--------------------

Use :py:class:`RCCLogicalVolume` for right-circular-cylinder regions.

Important constructor parameters:

* ``r``
* ``x0``, ``y0``, ``z0``: base point of the cylinder
* ``vx``, ``vy``, ``vz``: extrusion vector of the cylinder axis

This logical volume is useful for rods, pins, cylindrical detectors, beam-like
regions, and general cylinder-shaped selections that are not necessarily
aligned with one coordinate axis.

Examples:

.. code-block:: python

   # Cylinder aligned with z.
   pin = RCCLogicalVolume(
       r=0.4,
       x0=0.0, y0=0.0, z0=0.0,
       vx=0.0, vy=0.0, vz=1.5,
   )

   # Cylinder aligned with y.
   detector = RCCLogicalVolume(
       r=0.2,
       x0=0.0, y0=-0.5, z0=0.0,
       vx=0.0, vy=1.5, vz=0.0,
   )

Use this when a circular cross section plus an axis direction is the natural
description of the region.

``SphereLogicalVolume``
-----------------------

Use :py:class:`SphereLogicalVolume` for spherical regions.

Important constructor parameters:

* ``r``
* ``x``, ``y``, ``z``: sphere center

This is useful for spherical sources, spherical detectors, spherical exclusion
regions, or quick radial region definitions.

Example:

.. code-block:: python

   source_region = SphereLogicalVolume(
       r=0.5,
       x=0.0, y=0.0, z=0.0,
   )

``BooleanLogicalVolume``
------------------------

Use :py:class:`BooleanLogicalVolume` when a region is easiest to describe by
combining simpler logical volumes.

Its constructor takes a ``parts`` list. Each entry is a dictionary with:

* ``lv``: the logical volume to include in the boolean expression
* ``op``: ``True`` to include that part, ``False`` to exclude it

This is most useful for set-like constructions such as:

* a box with a cylindrical hole,
* a cylinder with a spherical cutout,
* "inside A but outside B."

Example:

.. code-block:: python

   outer = SphereLogicalVolume(r=1.0, x=0.0, y=0.0, z=0.0)
   hole = RCCLogicalVolume(
       r=0.2,
       x0=0.0, y0=0.0, z0=-1.0,
       vx=0.0, vy=0.0, vz=2.0,
   )
   region = BooleanLogicalVolume(
       parts=[
           {"op": True, "lv": outer},
           {"op": False, "lv": hole},
       ]
   )

Boolean logical volumes are especially valuable when the region of interest is
easy to reason about geometrically but awkward to describe with a single
primitive shape.

``SurfaceMeshLogicalVolume``
----------------------------

Use :py:class:`SurfaceMeshLogicalVolume` when the region should be defined by a
closed surface mesh rather than by one of the built-in primitive shapes.

Important constructor parameter:

* ``surface_mesh``: a :py:class:`pyopensn.mesh.SurfaceMesh`

This is the most flexible logical-volume type, but it is also the one that
depends most strongly on having a clean, well-defined input surface. It is
appropriate when the region boundary comes from a CAD-like surface description
and is not naturally a box, sphere, or cylinder.

Inspecting a logical volume directly
------------------------------------

All logical volumes derive from :py:class:`LogicalVolume` and expose:

* :py:meth:`LogicalVolume.Inside`

You can use this directly to check whether a point is inside a volume:

.. code-block:: python

   from pyopensn.math import Vector3

   lv = SphereLogicalVolume(r=0.5, x=0.0, y=0.0, z=0.0)
   print(lv.Inside(Vector3(0.1, 0.0, 0.0)))

This can be helpful when debugging a geometric selector before using it for
mesh labeling or source definition.

Practical Mesh Workflows
========================

Structured box-like geometry
----------------------------

Use:

* ``OrthogonalMeshGenerator``
* ``SetOrthogonalBoundaries``
* ``SetUniformBlockID`` or simple logical-volume block assignment

This is usually the easiest path for test problems and early input development.

Imported mesh
-------------

Use:

* ``FromFileMeshGenerator``
* a partitioner if needed
* a relabeling pass for blocks and boundaries

This is the normal path when geometry was created in another tool.

Extruded 3D transport mesh
--------------------------

Use:

* a base 2D mesh generator
* ``ExtruderMeshGenerator``
* then standard block and boundary labeling on the result

This is useful when the natural geometry description is layered.

Inspecting and Exporting Meshes
===============================

Use :py:meth:`MeshContinuum.ExportToPVTU` when you need to inspect the mesh,
the block-id layout, or the domain decomposition.

Practical Guidance
==================

As a practical guide:

* choose :py:class:`OrthogonalMeshGenerator` when the geometry is naturally
  box-like and easy to describe with node locations,
* choose :py:class:`FromFileMeshGenerator` when the mesh already exists in
  another tool,
* choose :py:class:`ExtruderMeshGenerator` when a 2D input geometry can 
  be cleanly extruded in the third dimension,
* use logical volumes to assign block ids, boundaries, or source regions after
  the mesh exists,
* start with a simple labeling scheme and only add complexity when the base
  problem is already behaving correctly.

.. note::

   A clean mesh workflow usually pays off more than a clever one. If a mesh,
   block-id map, and boundary-name set are all easy to explain, the rest of the
   transport input tends to be much easier to debug and maintain.
