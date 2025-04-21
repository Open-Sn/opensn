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
   :template: python

   math.Ylm

Point
^^^^^

.. autosummary::
   :toctree: generated
   :nosignatures:
   :template: python

   math.Vector3

Function wrappers
^^^^^^^^^^^^^^^^^

.. note::

   This API is temporary. It will no longer required once ``ParameterBlock`` is removed.

.. autosummary::
   :toctree: generated
   :nosignatures:
   :template: function

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
   :template: noinit

   aquad.AngularQuadrature

Product quadratures
^^^^^^^^^^^^^^^^^^^

.. autosummary::
   :toctree: generated
   :nosignatures:
   :template: noinit

   aquad.ProductQuadrature
   aquad.CurvilinearQuadrature

.. autosummary::
   :toctree: generated
   :nosignatures:
   :template: python

   aquad.GLProductQuadrature1DSlab
   aquad.GLCProductQuadrature2DXY
   aquad.GLCProductQuadrature3DXYZ

.. autosummary::
   :toctree: generated
   :nosignatures:
   :template: python

   aquad.GLCProductQuadrature2DRZ

Simplified LDFES quadrature
^^^^^^^^^^^^^^^^^^^^^^^^^^^

.. autosummary::
   :toctree: generated
   :nosignatures:
   :template: python

   aquad.SLDFESQuadrature


Field functions
---------------

Base class
^^^^^^^^^^

.. autosummary::
   :toctree: generated
   :nosignatures:
   :template: noinit

   fieldfunc.FieldFunction

Grid-based
^^^^^^^^^^

.. autosummary::
   :toctree: generated
   :nosignatures:
   :template: noinit

   fieldfunc.FieldFunctionGridBased

Interpolation
^^^^^^^^^^^^^

.. autosummary::
   :toctree: generated
   :nosignatures:
   :template: python

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
   :template: noinit

   mesh.MeshContinuum

Surface mesh
^^^^^^^^^^^^

.. autosummary::
   :toctree: generated
   :nosignatures:
   :template: python

   mesh.SurfaceMesh

Mesh generator
^^^^^^^^^^^^^^

.. autosummary::
   :toctree: generated
   :nosignatures:
   :template: noinit

   mesh.MeshGenerator

.. autosummary::
   :toctree: generated
   :nosignatures:
   :template: python

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
   :template: noinit

   mesh.GraphPartitioner

.. autosummary::
   :toctree: generated
   :nosignatures:
   :template: python

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
   :template: noinit

   logvol.LogicalVolume

Logical volume types
^^^^^^^^^^^^^^^^^^^^

.. autosummary::
   :toctree: generated
   :nosignatures:
   :template: python

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
   :template: python

   response.ResponseEvaluator


Source
------

.. autosummary::
   :toctree: generated
   :nosignatures:
   :template: python

   source.PointSource
   source.VolumetricSource


Cross section
-------------

.. autosummary::
   :toctree: generated
   :nosignatures:
   :template: python

   xs.MultiGroupXS


Problem
-------

Base class
^^^^^^^^^^

.. autosummary::
   :toctree: generated
   :nosignatures:
   :template: noinit

   solver.LBSProblem

Discrete ordinates problem
^^^^^^^^^^^^^^^^^^^^^^^^^^

.. autosummary::
   :toctree: generated
   :nosignatures:
   :template: python

   solver.DiscreteOrdinatesProblem

Discrete ordinates curvilinear problem
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

.. autosummary::
   :toctree: generated
   :nosignatures:
   :template: python.rst

   solver.DiscreteOrdinatesCurvilinearProblem


Solvers
-------

Solver base class
^^^^^^^^^^^^^^^^^

.. autosummary::
   :toctree: generated
   :nosignatures:
   :template: noinit.rst

   solver.Solver

Steady state solver
^^^^^^^^^^^^^^^^^^^

.. autosummary::
   :toctree: generated
   :nosignatures:
   :template: python

   solver.SteadyStateSolver

Non-linear k-eigen
^^^^^^^^^^^^^^^^^^

.. autosummary::
   :toctree: generated
   :nosignatures:
   :template: python

   solver.NonLinearKEigenSolver

Power iteration solver
^^^^^^^^^^^^^^^^^^^^^^

.. autosummary::
   :toctree: generated
   :nosignatures:
   :template: python

   solver.PowerIterationKEigenSolver
   solver.PowerIterationKEigenSCDSASolver
   solver.PowerIterationKEigenSMMSolver

Point kinetic transient solver
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

.. autosummary::
   :toctree: generated
   :nosignatures:
   :template: python

   solver.PRKSolver

Diffusion solver
^^^^^^^^^^^^^^^^

.. autosummary::
   :toctree: generated
   :nosignatures:
   :template: noinit

   diffusion.DiffusionSolverBase

.. autosummary::
   :toctree: generated
   :nosignatures:
   :template: python

   diffusion.CFEMDiffusionSolver
   diffusion.DFEMDiffusionSolver


Post-processors
---------------

Base class
^^^^^^^^^^

.. autosummary::
   :toctree: generated
   :nosignatures:
   :template: noinit

   post.PostProcessor

Post-processor
^^^^^^^^^^^^^^

.. autosummary::
   :toctree: generated
   :nosignatures:
   :template: python

   post.SolverInfoPostProcessor
   post.AggregateNodalValuePostProcessor
   post.CellVolumeIntegralPostProcessor

Printer
^^^^^^^

.. autosummary::
   :toctree: generated
   :nosignatures:
   :template: python

   post.Print
   post.SetPrinterOptions


Settings
--------

Logs
^^^^

.. autosummary::
   :toctree: generated
   :nosignatures:
   :template: noinit

   context.SetVerbosityLevel
   context.UseColor
   context.EnablePETScErrorHandler

Caliper configuration
^^^^^^^^^^^^^^^^^^^^^

.. autosummary::
   :toctree: generated
   :nosignatures:
   :template: noinit

   context.SetCaliperConfig
   context.EnableCaliper

Argument vector
^^^^^^^^^^^^^^^

.. autosummary::
   :toctree: generated
   :nosignatures:
   :template: noinit

   context.InitializeWithArgv
   context.Finalize
