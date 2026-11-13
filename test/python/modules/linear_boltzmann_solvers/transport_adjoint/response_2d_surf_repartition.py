#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
2D surface-response test with forward and adjoint data written on different mesh partitions.

This is the problem of response_2d_surf.py: the interior plane x = 6 separates the source
(x < 1) from the detector (x > 9.5) and all exterior boundaries are vacuum, so

    QoI = int_S int_4pi (Omega . n) psi(r, Omega) psi^dagger(r, Omega) dOmega dA,   n = +x,

equals the volumetric forward QoI and the adjoint response. For the upwind DFEM discretization
this identity is exact when each trace is taken from its upwind side.

The forward solution is computed on a mesh partitioned at x = 6, so the two sides of every face of
the plane are written by different ranks. The adjoint solution is computed on a mesh of the same
geometry partitioned in y. ResponseEvaluator.EvaluateSurfaceResponse must match the faces by
location across the two file sets and reproduce the volumetric QoIs. All three values are checked
against the response_2d_surf.py gold value.
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
    if size not in (2, 4):
        sys.exit(f"Expected 2 or 4 processors but got {size}.")

    num_cells = 20
    length = 10.0
    cell_size = length / num_cells
    nodes = [i * cell_size for i in range(num_cells + 1)]
    plane = 6.0

    # Forward partitions are cut at the interior plane; adjoint partitions are cut in y only.
    forward_xcuts = [plane] if size == 2 else [3.0, plane, 8.0]
    adjoint_ycuts = [5.0] if size == 2 else [2.5, 5.0, 7.5]

    source_volume = RPPLogicalVolume(
        xmin=nodes[0], xmax=nodes[2], ymin=nodes[0], ymax=nodes[2], infz=True
    )
    detector_volume = RPPLogicalVolume(
        xmin=nodes[-2], xmax=nodes[-1], ymin=nodes[-2], ymax=nodes[-1], infz=True
    )

    def make_grid(partitioner):
        grid = OrthogonalMeshGenerator(node_sets=[nodes, nodes], partitioner=partitioner).Execute()
        grid.SetOrthogonalBoundaries()
        grid.SetUniformBlockID(0)
        grid.SetBlockIDFromLogicalVolume(source_volume, 1, True)
        grid.SetBlockIDFromLogicalVolume(detector_volume, 2, True)
        return grid

    forward_grid = make_grid(
        KBAGraphPartitioner(nx=size, ny=1, nz=1, xcuts=forward_xcuts)
    )
    adjoint_grid = make_grid(
        KBAGraphPartitioner(nx=1, ny=size, nz=1, ycuts=adjoint_ycuts)
    )

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
    detector_sigma_a = 0.8
    source_strength = 1.0
    source_area = (nodes[2] - nodes[0]) ** 2

    quadrature = GLCProductQuadrature2DXY(n_polar=2, n_azimuthal=32, scattering_order=0)

    def make_problem(grid):
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

    prefix_base = f"Repartition{size}_"
    forward_surface_prefix = prefix_base + "FwdSurf_p"
    adjoint_surface_prefix = prefix_base + "AdjSurf_p"
    adjoint_flux_prefix = prefix_base + "AdjFluxMoments_p"
    interior_surfaces = {"inter_x": ("x", plane)}

    # Forward solve on the forward partition.
    forward_problem = make_problem(forward_grid)
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

    # Adjoint solve on the adjoint partition.
    adjoint_problem = make_problem(adjoint_grid)
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

    # The volumetric adjoint response reads the adjoint flux moments on their own partition. The
    # surface data sets come from two different partitions.
    evaluator = ResponseEvaluator(problem=adjoint_problem)
    # The interior surface name reads both stored sides, inter_x_u and inter_x_d.
    interior_tags = ["inter_x"]
    evaluator.SetOptions(
        sources={
            "material": [
                {"block_id": source_block_id, "strength": [source_strength / source_area]}
            ]
        },
        data=[
            {"name": "detector", "file_prefixes": {"flux_moments": adjoint_flux_prefix}},
            {
                "name": "forward",
                "file_prefixes": {"surface_angular_fluxes": forward_surface_prefix},
                "surfaces": interior_tags,
            },
            {
                "name": "adjoint",
                "file_prefixes": {"surface_angular_fluxes": adjoint_surface_prefix},
                "surfaces": interior_tags,
            },
        ],
    )
    adjoint_qoi = evaluator.EvaluateResponse("detector")
    (surface,) = evaluator.EvaluateSurfaceResponse("forward", "adjoint").values()
    surface_qoi = surface["total"]

    if rank == 0:
        print(f"Forward QoI={forward_qoi:.8e}")
        print(f"Adjoint QoI={adjoint_qoi:.8e}")
        print(f"Repartitioned Surface QoI={surface_qoi:.8e}")
        # Relative error of the surface identity: the surface QoI must equal the forward QoI.
        print(f"Repartitioned identity error={abs(surface_qoi / forward_qoi - 1.0):.3e}")

    barrier()
    for file_prefix in (forward_surface_prefix, adjoint_surface_prefix, adjoint_flux_prefix):
        os.remove(f"{file_prefix}{rank}.h5")
