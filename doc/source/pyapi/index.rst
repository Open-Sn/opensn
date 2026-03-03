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

Quadrature points
^^^^^^^^^^^^^^^^^

.. autosummary::
   :toctree: generated
   :nosignatures:
   :template: python.rst

   aquad.QuadraturePointPhiTheta

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

Triangular quadrature
^^^^^^^^^^^^^^^^^^^^^

.. autosummary::
   :toctree: generated
   :nosignatures:
   :template: noinit.rst

   aquad.TriangularQuadrature

.. autosummary::
   :toctree: generated
   :nosignatures:
   :template: python.rst

   aquad.GLCTriangularQuadrature2DXY
   aquad.GLCTriangularQuadrature3DXYZ

Lebedev quadrature
^^^^^^^^^^^^^^^^^^

.. autosummary::
   :toctree: generated
   :nosignatures:
   :template: python.rst

   aquad.LebedevQuadrature2DXY
   aquad.LebedevQuadrature3DXYZ

Simplified LDFES quadrature
^^^^^^^^^^^^^^^^^^^^^^^^^^^

.. autosummary::
   :toctree: generated
   :nosignatures:
   :template: python.rst

   aquad.SLDFEsqQuadrature3DXYZ
   aquad.SLDFEsqQuadrature2DXY


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

   fieldfunc.FieldFunctionInterpolation
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

.. note::

   **Forward/adjoint mode transitions are destructive by design.**
   In Python there are three valid ways to set mode:

   1. Set ``options={'adjoint': ...}`` in the problem constructor.
   2. Call ``problem.SetOptions(adjoint=...)``.
   3. Call ``problem.SetAdjoint(...)`` (low-level equivalent).

   A mode transition triggers:

   - reinitializes material mode (forward vs adjoint),
   - clears point and volumetric sources,
   - clears boundary conditions,
   - zeros scalar and angular flux state.

   The block-id to cross-section map is preserved. After switching mode, reapply
   the desired driving terms (sources and boundaries) before solving.

   ``SetOptions`` is additive: only explicitly supplied options are updated.
   If ``adjoint`` is omitted in a ``SetOptions`` call, the current mode is unchanged.

   Adjoint mode is applied to the mapped ``MultiGroupXS`` objects themselves.
   If the same cross-section object is shared across multiple problems, toggling
   adjoint mode in one problem affects all problems using that object.

Base class
^^^^^^^^^^

.. autosummary::
   :toctree: generated
   :nosignatures:
   :template: noinit.rst

   solver.Problem
   solver.LBSProblem

Discrete ordinates problem
^^^^^^^^^^^^^^^^^^^^^^^^^^

.. note::

   **Steady-state <-> transient mode transitions are non-destructive for problem setup.**
   Calling ``SetTimeDependentMode()`` or ``SetSteadyStateMode()`` preserves:

   - mesh and discretization,
   - block-id to cross-section mapping,
   - boundary conditions and source definitions,
   - scalar flux moments (``phi_new``/``phi_old`` state).

   For ``SetTimeDependentMode()``, OpenSn needs angular fluxes
   (``psi``) for the transient RHS time term. Transition behavior is:

   - If ``save_angular_flux`` is already enabled:
     - the mode switch performs internal transient mode configuration only.
   - If ``save_angular_flux`` is disabled:
     - OpenSn enables angular flux storage,
     - caches the scalar flux and performs a sweep to reconstruct ``psi``,
     - restores the cached scalar flux,
     - completes internal transient mode configuration.

   For ``SetSteadyStateMode()``, transient-only internals are reset. If angular
   flux saving had been temporarily forced for transient use, the original
   ``save_angular_flux`` user setting is restored.

   Practical implication:
   - ``save_angular_flux=False`` is valid for steady-state solves that later
     transition to transient; OpenSn will reconstruct transient ``psi`` as needed.

.. warning::

   Reconstructed angular flux (``psi``) is consistent with the converged scalar
   scalar flux moments, but it is not guaranteed to be identical to the angular
   flux that would have been available if ``save_angular_flux=True`` had been
   used during the steady-state solve. Small first-step transient differences
   can therefore appear, especially for problems with reflecting boundaries or
   lagged angular fluxes.

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

Transient solver
^^^^^^^^^^^^^^^^

.. autosummary::
   :toctree: generated
   :nosignatures:
   :template: python.rst

   solver.TransientSolver

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
