#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
2D transport response test exercising surface angular-flux I/O.
"""

import os
import sys

if "opensn_console" not in globals():
    from mpi4py import MPI

    size = MPI.COMM_WORLD.Get_size()
    rank = MPI.COMM_WORLD.Get_rank()
    barrier = MPI.COMM_WORLD.Barrier
    global_sum = MPI.COMM_WORLD.allreduce
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
    global_sum = MPIAllReduce


if __name__ == "__main__":
    num_procs = 1
    if size != num_procs:
        sys.exit(f"Incorrect number of processors. Expected {num_procs} processors but got {size}.")

    # Create mesh and assign material IDs.
    num_cells = 20
    length = 10.0
    cell_size = length / num_cells
    nodes = [i * cell_size for i in range(num_cells + 1)]
    meshgen = OrthogonalMeshGenerator(node_sets=[nodes, nodes])
    grid = meshgen.Execute()
    grid.SetOrthogonalBoundaries()
    grid.SetUniformBlockID(0)

    source_volume = RPPLogicalVolume(
        xmin=nodes[0], xmax=nodes[2], ymin=nodes[0], ymax=nodes[2], infz=True
    )
    grid.SetBlockIDFromLogicalVolume(source_volume, 1, True)

    detector_volume = RPPLogicalVolume(
        xmin=nodes[-2], xmax=nodes[-1], ymin=nodes[-2], ymax=nodes[-1], infz=True
    )
    grid.SetBlockIDFromLogicalVolume(detector_volume, 2, True)

    # Define cross sections.
    xs_background = MultiGroupXS()
    xs_background.CreateSimpleOneGroup(0.1, 0.9)

    xs_source = MultiGroupXS()
    xs_source.CreateSimpleOneGroup(2.0, 1.0)

    xs_detector = MultiGroupXS()
    xs_detector.CreateSimpleOneGroup(0.8, 0.0)

    xs_map = [
        {"block_ids": [0], "xs": xs_background},
        {"block_ids": [1], "xs": xs_source},
        {"block_ids": [2], "xs": xs_detector},
    ]
    source_block_id = 1

    # Set up the transport problem.
    quadrature = GLCProductQuadrature2DXY(n_polar=2, n_azimuthal=32, scattering_order=0)
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
                "gmres_restart_interval": 100,
            },
        ],
        xs_map=xs_map,
        boundary_conditions=[
            {"name": "xmin", "type": "vacuum"},
            {"name": "xmax", "type": "vacuum"},
            {"name": "ymin", "type": "vacuum"},
            {"name": "ymax", "type": "vacuum"},
        ],
        options={"save_angular_flux": True},
    )

    source_strength = 1.0
    source_area = (nodes[2] - nodes[0]) ** 2
    forward_source = VolumetricSource(
        block_ids=[source_block_id], group_strength=[source_strength / source_area]
    )
    problem.SetVolumetricSources(volumetric_sources=[forward_source])

    # Forward solve and surface-flux export.
    solver = SteadyStateSourceSolver(problem=problem)
    solver.Initialize()
    solver.Execute()

    field_functions = problem.GetScalarFluxFieldFunction()
    interpolator = FieldFunctionInterpolationVolume()
    interpolator.SetOperationType("sum")
    interpolator.SetLogicalVolume(detector_volume)
    interpolator.AddFieldFunction(field_functions[0])
    interpolator.Execute()
    forward_qoi = 0.8 * interpolator.GetValue()

    forward_surface_prefix = "InteriorSurf_FwdSrc_p"
    problem.WriteSurfaceAngularFluxes(forward_surface_prefix, [("inter_x", {"x": 6.0})])

    # Adjoint solve and response evaluation.
    problem.SetAdjoint(True)
    adjoint_source = VolumetricSource(logical_volume=detector_volume, group_strength=[0.8])
    problem.SetVolumetricSources(volumetric_sources=[adjoint_source])
    solver.Execute()

    adjoint_flux_prefix = "AdjFluxMoments_p"
    adjoint_surface_prefix = "InteriorSurf_AdjSrc_p"
    problem.WriteFluxMoments(adjoint_flux_prefix)
    problem.WriteSurfaceAngularFluxes(adjoint_surface_prefix, [("inter_x", {"x": 6.0})])

    evaluator = ResponseEvaluator(problem=problem)
    evaluator.SetOptions(
        buffers=[
            {"name": "detector", "file_prefixes": {"flux_moments": adjoint_flux_prefix}}
        ],
        sources={
            "material": [
                {"block_id": source_block_id, "strength": [source_strength / source_area]}
            ]
        },
    )
    adjoint_qoi = evaluator.EvaluateResponse("detector")

    forward_surface_up = problem.ReadSurfaceAngularFluxes(
        forward_surface_prefix, ["inter_x_u"]
    )
    forward_surface_down = problem.ReadSurfaceAngularFluxes(
        forward_surface_prefix, ["inter_x_d"]
    )
    local_forward_flux = sum(forward_surface_up[0]["data"]["psi"]) + sum(
        forward_surface_down[0]["data"]["psi"]
    )
    forward_boundary_flux = global_sum(local_forward_flux)

    adjoint_surface_up = problem.ReadSurfaceAngularFluxes(
        adjoint_surface_prefix, ["inter_x_u"]
    )
    adjoint_surface_down = problem.ReadSurfaceAngularFluxes(
        adjoint_surface_prefix, ["inter_x_d"]
    )
    local_adjoint_importance = sum(adjoint_surface_up[0]["data"]["psi"]) + sum(
        adjoint_surface_down[0]["data"]["psi"]
    )
    adjoint_boundary_importance = global_sum(local_adjoint_importance)

    if rank == 0:
        print(f"Forward QoI={forward_qoi:.5e}")
        print(f"Adjoint QoI={adjoint_qoi:.5e}")
        print(f"Forward Boundary Flux={forward_boundary_flux:.5e}")
        print(f"Adjoint Boundary Importance={adjoint_boundary_importance:.5e}")

    barrier()
    for file_prefix in (forward_surface_prefix, adjoint_flux_prefix, adjoint_surface_prefix):
        os.remove(f"{file_prefix}{rank}.h5")
