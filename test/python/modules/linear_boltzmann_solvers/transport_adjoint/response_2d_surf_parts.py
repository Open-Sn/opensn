#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
2D surface response parts with no return flux across the surface.

The interior plane x = 6 separates a pure-scattering source region (x < 1, 4 < y < 6) from a purely
absorbing detector facing it (x > 9.5, 4 < y < 6) in a void background, and all exterior boundaries
are vacuum. Nothing beyond the
plane emits or scatters, so no forward particle crosses the plane in the negative x direction: the
forward angular flux leaving the inter_x_d side is zero. The surface response across the plane,

    QoI = int_S int_4pi (Omega . n) psi(r, Omega) psi^dagger(r, Omega) dOmega dA,   n = +x,

then has only its positive-x part. This test checks that, for the interior surface inter_x,

* the total response equals the volumetric forward QoI;
* the outgoing part, carried in the positive x direction (the forward inter_x_u side against the
  adjoint inter_x_d side), also equals the volumetric forward QoI;
* the incoming part, carried in the negative x direction (the forward inter_x_d side against the
  adjoint inter_x_u side), is zero.

The mesh is partitioned at the plane, so the two sides of every face are on different ranks.
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
    from pyopensn.mesh import KBAGraphPartitioner, OrthogonalMeshGenerator
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

    num_cells = 20
    length = 10.0
    cell_size = length / num_cells
    nodes = [i * cell_size for i in range(num_cells + 1)]
    plane = 6.0

    # The source and detector face each other across the plane.
    source_volume = RPPLogicalVolume(
        xmin=nodes[0], xmax=nodes[2], ymin=nodes[8], ymax=nodes[12], infz=True
    )
    detector_volume = RPPLogicalVolume(
        xmin=nodes[-2], xmax=nodes[-1], ymin=nodes[8], ymax=nodes[12], infz=True
    )

    # Partition at the plane, so each side of the plane is on its own rank.
    grid = OrthogonalMeshGenerator(
        node_sets=[nodes, nodes],
        partitioner=KBAGraphPartitioner(nx=2, ny=1, nz=1, xcuts=[plane]),
    ).Execute()
    grid.SetOrthogonalBoundaries()
    grid.SetUniformBlockID(0)
    grid.SetBlockIDFromLogicalVolume(source_volume, 1, True)
    grid.SetBlockIDFromLogicalVolume(detector_volume, 2, True)

    # Void background, pure-scattering source region, purely absorbing detector.
    xs_void = MultiGroupXS()
    xs_void.CreateSimpleOneGroup(0.0, 0.0)
    xs_source = MultiGroupXS()
    xs_source.CreateSimpleOneGroup(2.0, 1.0)
    xs_detector = MultiGroupXS()
    xs_detector.CreateSimpleOneGroup(0.8, 0.0)
    xs_map = [
        {"block_ids": [0], "xs": xs_void},
        {"block_ids": [1], "xs": xs_source},
        {"block_ids": [2], "xs": xs_detector},
    ]
    source_block_id = 1
    detector_sigma_a = 0.8
    source_strength = 1.0
    source_area = (nodes[2] - nodes[0]) * (nodes[12] - nodes[8])

    quadrature = GLCProductQuadrature2DXY(n_polar=2, n_azimuthal=32, scattering_order=0)

    def make_problem():
        return DiscreteOrdinatesProblem(
            mesh=grid,
            num_groups=1,
            groupsets=[
                {
                    "groups_from_to": (0, 0),
                    "angular_quadrature": quadrature,
                    "inner_linear_method": "petsc_gmres",
                    "l_abs_tol": 1.0e-12,
                    "l_max_its": 500,
                    "gmres_restart_interval": 100,
                },
            ],
            xs_map=xs_map,
            boundary_conditions=[
                {"name": name, "type": "vacuum"} for name in ("xmin", "xmax", "ymin", "ymax")
            ],
            options={"save_angular_flux": True},
        )

    forward_surface_prefix = "SurfPartsFwdSurf_p"
    adjoint_surface_prefix = "SurfPartsAdjSurf_p"
    adjoint_flux_prefix = "SurfPartsAdjFluxMoments_p"
    interior_surfaces = {"inter_x": ("x", plane)}

    # Forward solve, volumetric QoI, and forward surface data.
    forward_problem = make_problem()
    forward_problem.SetVolumetricSources(
        volumetric_sources=[
            VolumetricSource(
                block_ids=[source_block_id], group_strength=[source_strength / source_area]
            )
        ]
    )
    forward_solver = SteadyStateSourceSolver(problem=forward_problem)
    forward_solver.Initialize()
    forward_solver.Execute()
    interpolator = FieldFunctionInterpolationVolume()
    interpolator.SetOperationType("sum")
    interpolator.SetLogicalVolume(detector_volume)
    interpolator.AddFieldFunction(forward_problem.GetScalarFluxFieldFunction()[0])
    interpolator.Execute()
    forward_qoi = detector_sigma_a * interpolator.GetValue()
    forward_problem.WriteSurfaceAngularFluxes(
        forward_surface_prefix, interior_surfaces=interior_surfaces
    )

    # Adjoint solve with the detector response as its source.
    adjoint_problem = make_problem()
    adjoint_problem.SetAdjoint(True)
    adjoint_problem.SetVolumetricSources(
        volumetric_sources=[
            VolumetricSource(logical_volume=detector_volume, group_strength=[detector_sigma_a])
        ]
    )
    adjoint_solver = SteadyStateSourceSolver(problem=adjoint_problem)
    adjoint_solver.Initialize()
    adjoint_solver.Execute()
    adjoint_problem.WriteFluxMoments(adjoint_flux_prefix)
    adjoint_problem.WriteSurfaceAngularFluxes(
        adjoint_surface_prefix, interior_surfaces=interior_surfaces
    )

    def surface_data(name, prefix, surface):
        return {
            "name": name,
            "file_prefixes": {"surface_angular_fluxes": prefix},
            "surfaces": [surface],
        }

    # "inter_x" reads both sides.
    evaluator = ResponseEvaluator(problem=adjoint_problem)
    evaluator.SetOptions(
        sources={
            "material": [
                {"block_id": source_block_id, "strength": [source_strength / source_area]}
            ]
        },
        data=[
            {"name": "detector", "file_prefixes": {"flux_moments": adjoint_flux_prefix}},
            surface_data("forward", forward_surface_prefix, "inter_x"),
            surface_data("adjoint", adjoint_surface_prefix, "inter_x"),
        ],
    )
    adjoint_qoi = evaluator.EvaluateResponse("detector")
    surface = evaluator.EvaluateSurfaceResponse("forward", "adjoint", ["inter_x"])["inter_x"]
    total_qoi = surface["total"]
    outgoing_qoi = surface["outgoing"]
    # The part carried in the negative x direction; incoming is its negative along +x.
    returning_qoi = -surface["incoming"]

    # A single side of an interior surface is rejected: a response needs both sides.
    lone_side_rejected = 0
    try:
        evaluator.SetOptions(data=[surface_data("lone", forward_surface_prefix, "inter_x_u")])
    except ValueError as error:
        lone_side_rejected = int("without its other side" in str(error))

    if rank == 0:
        print(f"Forward QoI={forward_qoi:.8e}")
        print(f"Adjoint QoI={adjoint_qoi:.8e}")
        print(f"Total QoI={total_qoi:.8e}")
        print(f"Outgoing part QoI={outgoing_qoi:.8e}")
        print(f"Returning part QoI={returning_qoi:.8e}")
        # Relative errors of the surface identity, and the returning part relative to the QoI.
        print(f"Total identity error={abs(total_qoi / forward_qoi - 1.0):.3e}")
        print(f"Outgoing part identity error={abs(outgoing_qoi / forward_qoi - 1.0):.3e}")
        print(f"Returning fraction={abs(returning_qoi / forward_qoi):.3e}")
        print(f"Lone side rejected={lone_side_rejected}")

    barrier()
    for file_prefix in (forward_surface_prefix, adjoint_surface_prefix, adjoint_flux_prefix):
        os.remove(f"{file_prefix}{rank}.h5")
