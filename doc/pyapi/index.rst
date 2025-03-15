Python API
==========

.. currentmodule:: pyopensn


Math
----

.. autosummary::
   :toctree: generated
   :nosignatures:
   :template: python.rst

   math.Ylm

.. autosummary::
   :toctree: generated
   :nosignatures:
   :template: python.rst

   math.Vector3

.. autosummary::
   :toctree: generated
   :nosignatures:
   :template: python.rst

   math.ScalarMaterialFunction
   math.ScalarSpatialMaterialFunction
   math.VectorSpatialFunction


Angular quadrature
------------------

.. autosummary::
   :toctree: generated
   :nosignatures:
   :template: python.rst

   aquad.AngularQuadrature
   aquad.ProductQuadrature
   aquad.GLProductQuadrature1DSlab
   aquad.GLCProductQuadrature2DXY
   aquad.GLCProductQuadrature3DXYZ
   aquad.CurvilinearQuadrature
   aquad.GLCProductQuadrature2DRZ
   aquad.SLDFESQuadrature


Field functions
---------------

.. autosummary::
   :toctree: generated
   :nosignatures:
   :template: python.rst

   fieldfunc.FieldFunction
   fieldfunc.FieldFunctionGridBased

.. autosummary::
   :toctree: generated
   :nosignatures:
   :template: python.rst

   fieldfunc.FieldFunctionInterpolation
   fieldfunc.FieldFunctionInterpolationPoint
   fieldfunc.FieldFunctionInterpolationLine
   fieldfunc.FieldFunctionInterpolationVolume


Mesh
----

.. autosummary::
   :toctree: generated
   :nosignatures:
   :template: python.rst

   mesh.MeshContinuum
   mesh.SurfaceMesh

.. autosummary::
   :toctree: generated
   :nosignatures:
   :template: python.rst

   mesh.MeshGenerator
   mesh.ExtruderMeshGenerator
   mesh.OrthogonalMeshGenerator
   mesh.FromFileMeshGenerator
   mesh.SplitFileMeshGenerator
   mesh.DistributedMeshGenerator

.. autosummary::
   :toctree: generated
   :nosignatures:
   :template: python.rst

   mesh.GraphPartitioner
   mesh.KBAGraphPartitioner
   mesh.LinearGraphPartitioner
   mesh.PETScGraphPartitioner


Logical volume
--------------

.. autosummary::
   :toctree: generated
   :nosignatures:
   :template: python.rst

   logvol.LogicalVolume
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

.. autosummary::
   :toctree: generated
   :nosignatures:
   :template: python.rst

   solver.Solver
   solver.LBSSolver
   solver.DiscreteOrdinatesSolver
   solver.DiscreteOrdinatesCurvilinearSolver
   solver.DiffusionDFEMSolver
   solver.SteadyStateSolver
   solver.NonLinearKEigen
   solver.PowerIterationKEigen
   solver.PowerIterationKEigenSCDSA
   solver.PowerIterationKEigenSMM
   solver.PRKSolver

.. autosummary::
   :toctree: generated
   :nosignatures:
   :template: python.rst

   diffusion.DiffusionSolverBase
   diffusion.CFEMDiffusionSolver
   diffusion.DFEMDiffusionSolver


Post-processors
---------------

.. autosummary::
   :toctree: generated
   :nosignatures:
   :template: python.rst

   post.PostProcessor
   post.SolverInfoPostProcessor
   post.AggregateNodalValuePostProcessor
   post.CellVolumeIntegralPostProcessor

.. autosummary::
   :toctree: generated
   :nosignatures:
   :template: python.rst

   post.Print
   post.SetPrinterOptions


Settings
--------

.. autosummary::
   :toctree: generated
   :nosignatures:
   :template: python.rst

   settings.SetVerbosityLevel
   settings.UseColor
   settings.EnablePETScErrorHandler
   settings.SetCaliperConfig
   settings.EnableCaliper

