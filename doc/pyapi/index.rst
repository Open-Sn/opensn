.. _pyapi:

Python API
==========

.. currentmodule:: pyopensn


Math
----

Spherical harmonics
^^^^^^^^^^^^^^^^^^^

.. autosummary::
   :toctree: generated
   :nosignatures:
   :template: python.rst

   math.Ylm

Point
^^^^^

.. autosummary::
   :toctree: generated
   :nosignatures:
   :template: python.rst

   math.Vector3

Function wrappers
^^^^^^^^^^^^^^^^^

.. note::

   This API is temporary. It will no longer required once ``ParameterBlock`` is removed.

.. autosummary::
   :toctree: generated
   :nosignatures:
   :template: function.rst

   math.Function
   math.ScalarMaterialFunction
   math.ScalarSpatialFunction
   math.ScalarSpatialMaterialFunction
   math.VectorSpatialFunction


Angular quadrature
------------------

Base class
^^^^^^^^^^

.. autosummary::
   :toctree: generated
   :nosignatures:
   :template: noinit.rst

   aquad.AngularQuadrature

Product quadratures
^^^^^^^^^^^^^^^^^^^

.. autosummary::
   :toctree: generated
   :nosignatures:
   :template: noinit.rst

   aquad.ProductQuadrature
   aquad.CurvilinearQuadrature

.. autosummary::
   :toctree: generated
   :nosignatures:
   :template: python.rst

   aquad.GLProductQuadrature1DSlab
   aquad.GLCProductQuadrature2DXY
   aquad.GLCProductQuadrature3DXYZ

.. autosummary::
   :toctree: generated
   :nosignatures:
   :template: python.rst

   aquad.GLCProductQuadrature2DRZ

Simplified LDFES quadrature
^^^^^^^^^^^^^^^^^^^^^^^^^^^

.. autosummary::
   :toctree: generated
   :nosignatures:
   :template: python.rst

   aquad.SLDFESQuadrature


Field functions
---------------

Base class
^^^^^^^^^^

.. autosummary::
   :toctree: generated
   :nosignatures:
   :template: noinit.rst

   fieldfunc.FieldFunction

Grid-based
^^^^^^^^^^

.. autosummary::
   :toctree: generated
   :nosignatures:
   :template: noinit.rst

   fieldfunc.FieldFunctionGridBased

Interpolation
^^^^^^^^^^^^^

.. autosummary::
   :toctree: generated
   :nosignatures:
   :template: python.rst

   fieldfunc.FieldFunctionInterpolationPoint
   fieldfunc.FieldFunctionInterpolationLine
   fieldfunc.FieldFunctionInterpolationVolume


Mesh
----

Mesh
^^^^

.. autosummary::
   :toctree: generated
   :nosignatures:
   :template: noinit.rst

   mesh.MeshContinuum

Surface mesh
^^^^^^^^^^^^

.. autosummary::
   :toctree: generated
   :nosignatures:
   :template: python.rst

   mesh.SurfaceMesh

Mesh generator
^^^^^^^^^^^^^^

.. autosummary::
   :toctree: generated
   :nosignatures:
   :template: noinit.rst

   mesh.MeshGenerator

.. autosummary::
   :toctree: generated
   :nosignatures:
   :template: python.rst

   mesh.ExtruderMeshGenerator
   mesh.OrthogonalMeshGenerator
   mesh.FromFileMeshGenerator
   mesh.SplitFileMeshGenerator
   mesh.DistributedMeshGenerator

Graph partitioner
^^^^^^^^^^^^^^^^^

.. autosummary::
   :toctree: generated
   :nosignatures:
   :template: noinit.rst

   mesh.GraphPartitioner

.. autosummary::
   :toctree: generated
   :nosignatures:
   :template: python.rst

   mesh.KBAGraphPartitioner
   mesh.LinearGraphPartitioner
   mesh.PETScGraphPartitioner


Logical volume
--------------

Base class
^^^^^^^^^^

.. autosummary::
   :toctree: generated
   :nosignatures:
   :template: noinit.rst

   logvol.LogicalVolume

Logical volume types
^^^^^^^^^^^^^^^^^^^^

.. autosummary::
   :toctree: generated
   :nosignatures:
   :template: python.rst

   logvol.BooleanLogicalVolume
   logvol.RCCLogicalVolume
   logvol.RPPLogicalVolume
   logvol.SphereLogicalVolume
   logvol.SurfaceMeshLogicalVolume


Response evaluator
------------------

.. autosummary::
   :toctree: generated
   :nosignatures:
   :template: python.rst

   response.ResponseEvaluator


Source
------

.. autosummary::
   :toctree: generated
   :nosignatures:
   :template: python.rst

   source.PointSource
   source.VolumetricSource


Cross section
-------------

.. autosummary::
   :toctree: generated
   :nosignatures:
   :template: python.rst

   xs.MultiGroupXS


Solver
------

Base class
^^^^^^^^^^

.. autosummary::
   :toctree: generated
   :nosignatures:
   :template: noinit.rst

   solver.Solver

Linear Boltzmann solver
^^^^^^^^^^^^^^^^^^^^^^^

.. autosummary::
   :toctree: generated
   :nosignatures:
   :template: noinit.rst

   solver.LBSProblem

.. autosummary::
   :toctree: generated
   :nosignatures:
   :template: python.rst

   solver.DiscreteOrdinatesProblem
   solver.DiscreteOrdinatesCurvilinearProblem
   solver.DiffusionDFEMSolver

Steady state solver
^^^^^^^^^^^^^^^^^^^

.. autosummary::
   :toctree: generated
   :nosignatures:
   :template: python.rst

   solver.SteadyStateSolver

Non-linear k-eigen
^^^^^^^^^^^^^^^^^^

.. autosummary::
   :toctree: generated
   :nosignatures:
   :template: python.rst

   solver.NonLinearKEigenSolver

Power iteration solver
^^^^^^^^^^^^^^^^^^^^^^

.. autosummary::
   :toctree: generated
   :nosignatures:
   :template: python.rst

   solver.PowerIterationKEigenSolver
   solver.PowerIterationKEigenSCDSASolver
   solver.PowerIterationKEigenSMM

Point kinetic trnasient solver
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

.. autosummary::
   :toctree: generated
   :nosignatures:
   :template: python.rst

   solver.PRKSolver

Diffusion solver
^^^^^^^^^^^^^^^^

.. autosummary::
   :toctree: generated
   :nosignatures:
   :template: noinit.rst

   diffusion.DiffusionSolverBase

.. autosummary::
   :toctree: generated
   :nosignatures:
   :template: python.rst

   diffusion.CFEMDiffusionSolver
   diffusion.DFEMDiffusionSolver


Post-processors
---------------

Base class
^^^^^^^^^^

.. autosummary::
   :toctree: generated
   :nosignatures:
   :template: noinit.rst

   post.PostProcessor

Post-processor
^^^^^^^^^^^^^^

.. autosummary::
   :toctree: generated
   :nosignatures:
   :template: python.rst

   post.SolverInfoPostProcessor
   post.AggregateNodalValuePostProcessor
   post.CellVolumeIntegralPostProcessor

Printer
^^^^^^^

.. autosummary::
   :toctree: generated
   :nosignatures:
   :template: python.rst

   post.Print
   post.SetPrinterOptions


Settings
--------

Logs
^^^^

.. autosummary::
   :toctree: generated
   :nosignatures:
   :template: noinit.rst

   settings.SetVerbosityLevel
   settings.UseColor
   settings.EnablePETScErrorHandler

Caliper configuration
^^^^^^^^^^^^^^^^^^^^^

.. autosummary::
   :toctree: generated
   :nosignatures:
   :template: noinit.rst

   settings.SetCaliperConfig
   settings.EnableCaliper
