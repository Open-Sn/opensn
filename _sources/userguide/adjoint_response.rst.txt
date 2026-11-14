=========================================
Adjoint Problems and Response Evaluation
=========================================

An adjoint solution gives a response, such as a detector reaction rate, for any
forward source with one inner product instead of one forward solve per source.
It is also the input to the cross-section sensitivities of
:ref:`userguide_xs_sensitivity`. The equations and the discrete conventions are
in :doc:`../theory/adjoint`.

Supported Cases
===============

Adjoint mode is available for :py:class:`pyopensn.solver.DiscreteOrdinatesProblem`
with:

* Cartesian geometry (not :py:class:`pyopensn.solver.DiscreteOrdinatesCurvilinearProblem`),
* standard angular operators (``operator_method="standard"``) and quadratures
  that contain, for every direction, an opposite direction of equal weight,
* steady-state mode, solved with
  :py:class:`pyopensn.solver.SteadyStateSourceSolver`,
  :py:class:`pyopensn.solver.PowerIterationKEigenSolver`, or
  :py:class:`pyopensn.solver.NonLinearKEigenSolver`.

Time-dependent mode, CSDA (``csda_enabled``), and an uncollided flux are rejected
in adjoint mode. Groupsets, DSA, CMFD, and delayed neutrons
(``use_precursors``) work as in forward mode. The energy transfer and fission
production, including the delayed production, are transposed. The diffusion
coefficients of DSA, CMFD, and SMM use the transport correction of the
transposed P1 scattering. Precursor concentrations are not computed from an
adjoint flux and stay zero.

Selecting the Mode
==================

Select the mode with ``options={"adjoint": True}`` when the problem is built,
or switch an existing problem with :py:meth:`pyopensn.solver.LBSProblem.SetAdjoint`
(:py:meth:`pyopensn.solver.LBSProblem.SetForward` switches back).

.. warning::

   Changing the mode with ``SetAdjoint`` or ``SetForward`` clears **all point
   and volumetric sources**, **all source moments**, and **all boundary
   sources** (isotropic and arbitrary boundaries, which become vacuum), and
   zeroes the flux. A forward source is not an adjoint source, so apply the
   sources for the new mode before solving with
   :py:meth:`pyopensn.solver.LBSProblem.SetPointSources`,
   :py:meth:`pyopensn.solver.LBSProblem.SetVolumetricSources`,
   :py:meth:`pyopensn.solver.DiscreteOrdinatesProblem.SetBoundaryOptions`, or,
   with ``use_source_moments=True``,
   :py:meth:`pyopensn.solver.LBSProblem.ReadSourceMoments`. Otherwise the solve
   uses no sources, or raises an error when it needs source moments. Vacuum and reflecting boundaries are
   kept because they mean the same in both modes. OpenSn logs a warning on
   every mode change.

.. code-block:: python

   problem.SetAdjoint(True)
   problem.SetVolumetricSources(volumetric_sources=[detector_source])

Adjoint Sources and Boundary Conditions
=======================================

The adjoint source is the response function. For a reaction rate
:math:`R=\int\sum_g\sigma^g_{\mathrm{det}}\phi^g\,dV` over a detector region,
use a volumetric source with ``group_strength`` equal to
:math:`\sigma_{\mathrm{det}}` on that region. A point source with ``strength``
equal to :math:`\sigma_{\mathrm{det}}` gives the response at that point
instead. As in the forward problem, group values are group averages (see
:doc:`../theory/adjoint`).

Vacuum and reflecting boundaries mean the same in both modes. An isotropic
boundary source on an adjoint problem gives the response
:math:`\sum_{n:\,\Omega_n\cdot n>0} w_n\,\Omega_n\cdot n\,f^g\,\psi^g_n`
integrated over that boundary, the outgoing partial current weighted by the
group strength :math:`f^g`.

OpenSn solves the adjoint with the forward sweeps and transposed cross sections.
Direction-dependent adjoint data, such as an ``"arbitrary"`` boundary callback
or source moments of odd degree, are therefore evaluated at the sweep direction
:math:`\Omega_n` but describe the adjoint for the opposite direction
:math:`-\Omega_n` (for 1D and 2D quadratures, the paired opposite direction
described in :doc:`../theory/adjoint`). Supply them accordingly. Isotropic data
need no change. Source moments made from an adjoint flux by
:py:meth:`pyopensn.solver.LBSProblem.CreateAndWriteSourceMoments` or
:py:meth:`pyopensn.solver.LBSProblem.ReadFluxMomentsAndMakeSourceMoments` are
already in this orientation.

Solving and Writing the Adjoint
===============================

The steady-state and k-eigenvalue solvers reorient the adjoint after the solve,
so the stored angular flux at direction index :math:`n` is the physical adjoint
for :math:`\Omega_n` and odd-degree flux moments have the physical sign.
:py:meth:`pyopensn.solver.LBSProblem.WriteFluxMoments` and
:py:meth:`pyopensn.solver.DiscreteOrdinatesProblem.WriteAngularFluxes` record
whether a file holds a forward or an adjoint solution and, for an adjoint, its
orientation. Write adjoint files after the solver has finished. A file written
before then is in solver orientation and is rejected by the response evaluator
and the sensitivity postprocessor.

Response Evaluation
===================

:py:class:`pyopensn.response.ResponseEvaluator` folds forward sources with
stored solutions, called data sets:

.. math::

   R = \int \phi^\dagger q\,dV
     + W \sum_{n:\,\Omega_n\cdot n < 0} w_n \int_{\partial D}
       |\Omega_n\cdot n|\,\psi^\dagger_n\,\psi^{\mathrm{in}}_n\,dA ,

over the material, point, volumetric, and isotropic boundary sources defined in
the evaluator. The boundary sum runs over the incoming directions of each
boundary face, :math:`w_n` is the quadrature weight, and :math:`W` is the sum of
the quadrature weights (1 for OpenSn quadratures, applied per groupset).

For an adjoint data set from a fixed-source solve
(:py:class:`pyopensn.solver.SteadyStateSourceSolver`), this equals the forward
response to the adjoint source that produced it, to the solver tolerance. The
forward problem must have the same mesh, cross sections, quadrature, and vacuum
and reflecting boundaries as the adjoint problem. A k-eigenvalue adjoint is the
fundamental-mode importance with an arbitrary normalization, not the solution
for a detector source, so folding it with a forward source gives no detector
response. It is used for k-eigenvalue sensitivities. The evaluator applies the
sources as the forward solve does:

* Volumetric sources are integrated with the consistent mass matrix, as the
  forward solve applies them.
* Point and volumetric sources are evaluated at the problem's time, and sources
  inactive then are skipped, as in a forward solve.
* Boundary sources need adjoint angular fluxes, so the adjoint problem needs
  ``save_angular_flux=True``.
* In RZ geometry the result includes the :math:`2\pi` azimuthal factor.

Flux-moment and angular-flux files must have been written with the same rank
count and partition as the problem, and angular-flux files with its quadrature.
A data set may also hold a forward solution, for example to integrate a forward
flux against a material source.

The evaluator also pairs forward and adjoint angular fluxes on boundaries and
interior surfaces written with
:py:meth:`pyopensn.solver.DiscreteOrdinatesProblem.WriteSurfaceAngularFluxes`.
These data sets use the ``surface_angular_fluxes`` prefix and are evaluated with
:py:meth:`pyopensn.response.ResponseEvaluator.EvaluateSurfaceResponse`, which
gives the response carried across each surface. Surface files may come from any
rank count and partition. See
:doc:`../tutorials/applications/detector_response/surface_response/surface_resp`.

Example
-------

.. code-block:: python

   from pyopensn.aquad import GLCProductQuadrature2DXY
   from pyopensn.logvol import RPPLogicalVolume
   from pyopensn.mesh import OrthogonalMeshGenerator
   from pyopensn.response import ResponseEvaluator
   from pyopensn.solver import DiscreteOrdinatesProblem, SteadyStateSourceSolver
   from pyopensn.source import VolumetricSource
   from pyopensn.xs import MultiGroupXS

   nodes = [i / 8 for i in range(9)]
   grid = OrthogonalMeshGenerator(node_sets=[nodes, nodes]).Execute()
   grid.SetUniformBlockID(0)
   xs = MultiGroupXS()
   xs.CreateSimpleOneGroup(1.0, 0.5)
   detector = RPPLogicalVolume(xmin=0.5, xmax=1.0, ymin=0.5, ymax=1.0, infz=True)

   # Adjoint problem: the source is the detector response function.
   problem = DiscreteOrdinatesProblem(
       mesh=grid,
       num_groups=1,
       groupsets=[
           {
               "groups_from_to": (0, 0),
               "angular_quadrature": GLCProductQuadrature2DXY(
                   n_polar=4, n_azimuthal=8, scattering_order=0
               ),
               "inner_linear_method": "petsc_gmres",
               "l_abs_tol": 1.0e-10,
           }
       ],
       xs_map=[{"block_ids": [0], "xs": xs}],
       volumetric_sources=[VolumetricSource(logical_volume=detector, group_strength=[1.0])],
       boundary_conditions=[{"name": "xmin", "type": "reflecting"}],
       options={"adjoint": True, "save_angular_flux": True},
   )
   solver = SteadyStateSourceSolver(problem=problem)
   solver.Initialize()
   solver.Execute()
   problem.WriteFluxMoments("detector_phi")
   problem.WriteAngularFluxes("detector_psi")

   # Responses for several forward source configurations.
   evaluator = ResponseEvaluator(problem=problem)
   evaluator.SetOptions(
       data=[
           {
               "name": "detector",
               "file_prefixes": {"flux_moments": "detector_phi", "angular_fluxes": "detector_psi"},
           }
       ]
   )
   for strength in (1.0, 2.0):
       evaluator.ClearForwardSources()
       evaluator.SetOptions(
           sources={
               "material": [{"block_id": 0, "strength": [strength]}],
               "boundary": [{"name": "ymin", "type": "isotropic", "group_strength": [0.5]}],
           }
       )
       response = evaluator.EvaluateResponse("detector")

Each response equals the detector reaction rate of a forward solve with a
uniform source of the given strength, an isotropic incident flux of 0.5 on
``ymin``, and a reflecting ``xmin``. ``EvaluateResponse`` is collective and
returns the same value on every rank.
