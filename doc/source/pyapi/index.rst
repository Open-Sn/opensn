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

.. autosummary::
   :toctree: generated
   :nosignatures:
   :template: function.rst

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
   aquad.CurvilinearProductQuadrature

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


Problem
-------

Base class
^^^^^^^^^^

.. autosummary::
   :toctree: generated
   :nosignatures:
   :template: noinit.rst

   solver.LBSProblem

Discrete ordinates problem
^^^^^^^^^^^^^^^^^^^^^^^^^^

.. autosummary::
   :toctree: generated
   :nosignatures:
   :template: python.rst

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

Steady state source solver
^^^^^^^^^^^^^^^^^^^^^^^^^^

.. autosummary::
   :toctree: generated
   :nosignatures:
   :template: python.rst

   solver.SteadyStateSourceSolver

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


Acceleration methods
--------------------

Discrete ordinates k-eigen acceleration base class
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

.. autosummary::
   :toctree: generated
   :nosignatures:
   :template: noinit.rst

   solver.DiscreteOrdinatesKEigenAcceleration

Discrete ordinates k-eigen acceleration
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

.. autosummary::
   :toctree: generated
   :nosignatures:
   :template: python.rst

   solver.SCDSAAcceleration
   solver.SMMAcceleration


Settings
--------

.. important::

   Functions in this section are only available to the module mode. For the
   console mode, refer to ``opensn --help`` for more information.

Logs
^^^^

.. autosummary::
   :toctree: generated
   :nosignatures:
   :template: noinit.rst

   context.SetVerbosityLevel
   context.UseColor
   context.EnablePETScErrorHandler

Caliper configuration
^^^^^^^^^^^^^^^^^^^^^

.. autosummary::
   :toctree: generated
   :nosignatures:
   :template: noinit.rst

   context.SetCaliperConfig
   context.EnableCaliper

Argument vector
^^^^^^^^^^^^^^^

.. autosummary::
   :toctree: generated
   :nosignatures:
   :template: noinit.rst

   context.InitializeWithArgv
   context.Finalize

GPU configuration
^^^^^^^^^^^^^^^^^

.. important::

   This API is only available when ``pyopensn.can_support_gpus`` is true.

.. note::

   The GPU assignment functionality is intended for use on single-user workstations.

   It should not be used for production runs on large systems. When multiple GPUs are available on a
   node, the mapping between MPI ranks and GPUs should be handled by the job scheduler.

.. autosummary::
   :toctree: generated
   :nosignatures:
   :template: noinit.rst

   device.get_device_count
   device.get_current_device
   device.set_device
