Introduction to Meshing and Partitioning
========================================

Meshing
-------

Python mesh construction starts with :py:class:`~pyopensn.mesh.MeshGenerator`. Calling
:py:meth:`~pyopensn.mesh.MeshGenerator.Execute` completes the generator pipeline and returns the
:py:class:`~pyopensn.mesh.MeshContinuum` used by a problem.

The ``pyopensn.mesh`` module provides these mesh generators:

* :py:class:`~pyopensn.mesh.OrthogonalMeshGenerator` creates one-, two-, or three-dimensional
  orthogonal meshes from coordinate-node sets.
* :py:class:`~pyopensn.mesh.FromFileMeshGenerator` imports an externally generated mesh.
* :py:class:`~pyopensn.mesh.ExtruderMeshGenerator` turns a two-dimensional input mesh into a
  layered three-dimensional mesh.
* :py:class:`~pyopensn.mesh.SplitFileMeshGenerator` writes or reads rank-specific binary mesh
  partitions for reuse.
* :py:class:`~pyopensn.mesh.DistributedMeshGenerator` partitions a mesh on rank zero and
  distributes the resulting partitions to the MPI ranks.

:py:class:`~pyopensn.mesh.FromFileMeshGenerator` accepts meshes in the following formats:

* ``msh`` (Gmsh)
* ``vtu`` and ``pvtu`` (VTK unstructured grids)
* ``e`` (Exodus II)
* ``case`` (Ensight Gold)
* ``obj`` (Wavefront surface meshes, for two-dimensional meshing)

Imported meshes can carry block IDs for materials and volumetric sources and boundary IDs for
boundary conditions and sources. IDs can also be assigned after generation using the methods on
:py:class:`~pyopensn.mesh.MeshContinuum`.

Mesh generators can be chained through their ``inputs`` parameter, allowing one generator to
transform the output of another.

Partitioning
------------

For parallel simulations, each mesh generator accepts a
:py:class:`~pyopensn.mesh.GraphPartitioner`. OpenSn provides:

* :py:class:`~pyopensn.mesh.PETScGraphPartitioner`, using either ParMETIS or PT-Scotch through
  PETSc. This is the default, with ParMETIS selected unless configured otherwise.
* :py:class:`~pyopensn.mesh.KBAGraphPartitioner`, which places cuts along coordinate axes and is
  especially useful for orthogonal grids.
* :py:class:`~pyopensn.mesh.LinearGraphPartitioner`, a simple index-based partitioner intended
  primarily for tests.

See the :doc:`Python API reference </pyapi/index>` for all generator and partitioner parameters.
