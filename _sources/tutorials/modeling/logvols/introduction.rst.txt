Introduction to Logical Volumes
===============================

A logical volume is a geometric selector. Its
:py:meth:`~pyopensn.logvol.LogicalVolume.Inside` method determines whether a point lies inside
the selected region. OpenSn can use logical volumes to assign mesh block or boundary IDs and to
identify subregions for sources and post-processing.

The ``pyopensn.logvol`` module provides the following concrete logical volumes:

* :py:class:`~pyopensn.logvol.RPPLogicalVolume` selects an axis-aligned rectangular
  parallelepiped. Individual axes can be made infinite, which also makes this class useful for
  slabs and extruded rectangles.
* :py:class:`~pyopensn.logvol.RCCLogicalVolume` selects a finite right circular cylinder defined
  by a base point, an axis vector, and a radius.
* :py:class:`~pyopensn.logvol.SphereLogicalVolume` selects a sphere defined by its center and
  radius.
* :py:class:`~pyopensn.logvol.BooleanLogicalVolume` combines other logical volumes with inclusion
  and exclusion operations.
* :py:class:`~pyopensn.logvol.SurfaceMeshLogicalVolume` selects the region enclosed by a
  :py:class:`~pyopensn.mesh.SurfaceMesh`.

All of these derive from :py:class:`~pyopensn.logvol.LogicalVolume`. See the
:doc:`Python API reference </pyapi/index>` for their constructor parameters and methods.

When a logical volume is passed to
:py:meth:`~pyopensn.mesh.MeshContinuum.SetBlockIDFromLogicalVolume`, OpenSn selects a cell by
testing its centroid. Boundary IDs can similarly be assigned with
:py:meth:`~pyopensn.mesh.MeshContinuum.SetBoundaryIDFromLogicalVolume`. A logical volume that
only intersects a cell, without containing its centroid, therefore does not select that cell.
