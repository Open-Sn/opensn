#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
2D adjoint response to an isotropic boundary source.

By forward-adjoint duality, the detector response to an isotropic incoming flux on xmin,

    QoI = int_V sigma_d phi dV,

equals the adjoint evaluation

    QoI = int_{xmin} int_{Omega . n < 0} |Omega . n| psi_b psi^dagger(r, Omega) dOmega dA.

The adjoint evaluation reads the stored adjoint angular flux at incoming directions, so it
requires the adjoint solution to be reoriented to the physical direction. 2D quadratures contain
only directions with Omega_z > 0, which makes this the case that exercises the in-plane direction
pairing.
"""

import os
import sys

if "opensn_console" not in globals():
    from mpi4py import MPI

    size = MPI.COMM_WORLD.Get_size()
    rank = MPI.COMM_WORLD.Get_rank()
    barrier = MPI.COMM_WORLD.Barrier
    sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), "../../../../../")))
    from pyopensn.aquad import GLCProductQuadrature2DXY
    from pyopensn.fieldfunc import FieldFunctionInterpolationVolume
    from pyopensn.logvol import RPPLogicalVolume
    from pyopensn.mesh import OrthogonalMeshGenerator
    from pyopensn.response import ResponseEvaluator
    from pyopensn.solver import DiscreteOrdinatesProblem, SteadyStateSourceSolver
    from pyopensn.source import VolumetricSource
    from pyopensn.xs import MultiGroupXS
else:
    barrier = MPIBarrier


if __name__ == "__main__":
    num_procs = 2
    if size != num_procs:
        sys.exit(f"Incorrect number of processors. Expected {num_procs} processors but got {size}.")

    nodes = [i * 0.5 for i in range(21)]
    grid = OrthogonalMeshGenerator(node_sets=[nodes, nodes]).Execute()
    grid.SetOrthogonalBoundaries()
    grid.SetUniformBlockID(0)

    detector_volume = RPPLogicalVolume(xmin=9.0, xmax=10.0, ymin=4.0, ymax=6.0, infz=True)
    grid.SetBlockIDFromLogicalVolume(detector_volume, 1, True)

    xs_background = MultiGroupXS()
    xs_background.CreateSimpleOneGroup(0.3, 0.5)
    xs_detector = MultiGroupXS()
    xs_detector.CreateSimpleOneGroup(0.8, 0.0)
    detector_sigma = 0.8
    boundary_strength = 1.0

    quadrature = GLCProductQuadrature2DXY(n_polar=4, n_azimuthal=16, scattering_order=0)

    def boundary_conditions(with_source):
        bcs = [{"name": name, "type": "vacuum"} for name in ("xmax", "ymin", "ymax")]
        if with_source:
            bcs.append(
                {"name": "xmin", "type": "isotropic", "group_strength": [boundary_strength]}
            )
        else:
            bcs.append({"name": "xmin", "type": "vacuum"})
        return bcs

    problem = DiscreteOrdinatesProblem(
        mesh=grid,
        num_groups=1,
        groupsets=[
            {
                "groups_from_to": (0, 0),
                "angular_quadrature": quadrature,
                "inner_linear_method": "petsc_gmres",
                "l_abs_tol": 1.0e-12,
                "l_max_its": 500,
            }
        ],
        xs_map=[
            {"block_ids": [0], "xs": xs_background},
            {"block_ids": [1], "xs": xs_detector},
        ],
        boundary_conditions=boundary_conditions(True),
        options={"save_angular_flux": True},
    )

    # Forward solve and volumetric detector response.
    solver = SteadyStateSourceSolver(problem=problem)
    solver.Initialize()
    solver.Execute()

    interpolator = FieldFunctionInterpolationVolume()
    interpolator.SetOperationType("sum")
    interpolator.SetLogicalVolume(detector_volume)
    interpolator.AddFieldFunction(problem.GetScalarFluxFieldFunction()[0])
    interpolator.Execute()
    forward_qoi = detector_sigma * interpolator.GetValue()

    # Adjoint solve and boundary-source response evaluation.
    problem.SetAdjoint(True)
    problem.SetBoundaryOptions(boundary_conditions=boundary_conditions(False))
    problem.SetVolumetricSources(
        volumetric_sources=[
            VolumetricSource(logical_volume=detector_volume, group_strength=[detector_sigma])
        ]
    )
    solver.Execute()

    adjoint_phi_prefix = "BndrySrcAdjPhi_p"
    adjoint_psi_prefix = "BndrySrcAdjPsi_p"
    problem.WriteFluxMoments(adjoint_phi_prefix)
    problem.WriteAngularFluxes(adjoint_psi_prefix)

    evaluator = ResponseEvaluator(problem=problem)
    evaluator.SetOptions(
        buffers=[
            {
                "name": "detector",
                "file_prefixes": {
                    "flux_moments": adjoint_phi_prefix,
                    "angular_fluxes": adjoint_psi_prefix,
                },
            }
        ],
        sources={
            "boundary": [
                {"name": "xmin", "type": "isotropic", "group_strength": [boundary_strength]}
            ]
        },
    )
    adjoint_qoi = evaluator.EvaluateResponse("detector")

    if rank == 0:
        print(f"Forward QoI={forward_qoi:.8e}")
        print(f"Adjoint QoI={adjoint_qoi:.8e}")

    barrier()
    for prefix in (adjoint_phi_prefix, adjoint_psi_prefix):
        os.remove(f"{prefix}{rank}.h5")
