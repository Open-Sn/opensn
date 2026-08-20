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
    num_procs = 2
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

    def make_problem(save_angular_flux):
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
                {"name": "xmin", "type": "vacuum"},
                {"name": "xmax", "type": "vacuum"},
                {"name": "ymin", "type": "vacuum"},
                {"name": "ymax", "type": "vacuum"},
            ],
            options={"save_angular_flux": save_angular_flux},
        )

    no_angular_flux_problem = make_problem(False)
    missing_storage_prefix = "MissingAngularFluxStorage_p"
    missing_storage_rejected = False
    try:
        no_angular_flux_problem.WriteSurfaceAngularFluxes(
            missing_storage_prefix, boundary_surfaces=["xmin"]
        )
    except ValueError as error:
        missing_storage_rejected = "options.save_angular_flux=true" in str(error)
    if not missing_storage_rejected:
        raise RuntimeError("Surface export without saved angular flux was not rejected")
    if os.path.exists(f"{missing_storage_prefix}{rank}.h5"):
        raise RuntimeError("Rejected surface export created an output file")

    problem = make_problem(True)

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

    invalid_boundary_prefix = "InvalidBoundary_p"
    invalid_boundary_file = f"{invalid_boundary_prefix}{rank}.h5"
    preserved_contents = b"preserve existing surface output"
    with open(invalid_boundary_file, "wb") as output:
        output.write(preserved_contents)
    try:
        invalid_boundary_rejected = False
        invalid_boundaries = ["missing"] if rank == 0 else []
        expected_error = (
            "not found in the boundary-name map"
            if rank == 0
            else "requested on another rank"
        )
        try:
            problem.WriteSurfaceAngularFluxes(
                invalid_boundary_prefix, boundary_surfaces=invalid_boundaries
            )
        except ValueError as error:
            invalid_boundary_rejected = expected_error in str(error)
        if not invalid_boundary_rejected:
            raise RuntimeError("Invalid boundary surface was not rejected")
        with open(invalid_boundary_file, "rb") as output:
            if output.read() != preserved_contents:
                raise RuntimeError("Rejected surface export modified an existing output file")
    finally:
        if os.path.exists(invalid_boundary_file):
            os.remove(invalid_boundary_file)

    field_functions = problem.GetScalarFluxFieldFunction()
    interpolator = FieldFunctionInterpolationVolume()
    interpolator.SetOperationType("sum")
    interpolator.SetLogicalVolume(detector_volume)
    interpolator.AddFieldFunction(field_functions[0])
    interpolator.Execute()
    forward_qoi = 0.8 * interpolator.GetValue()

    forward_surface_prefix = "InteriorSurf_FwdSrc_p"
    invalid_surface_prefix = "BoundaryAsInterior_p"
    boundary_as_internal_rejected = False
    try:
        problem.WriteSurfaceAngularFluxes(
            invalid_surface_prefix,
            boundary_surfaces=["xmin"],
            interior_surfaces={"not_interior": ("x", 0.0)},
        )
    except ValueError as error:
        boundary_as_internal_rejected = "coincides with an exterior boundary" in str(error)
    if not boundary_as_internal_rejected:
        raise RuntimeError("An interior surface coincident with a boundary was not rejected")
    if os.path.exists(f"{invalid_surface_prefix}{rank}.h5"):
        raise RuntimeError("Rejected surface selection created an output file")

    problem.WriteSurfaceAngularFluxes(
        forward_surface_prefix,
        boundary_surfaces=["xmin", "xmax", "ymin", "ymax"],
        interior_surfaces={"inter_x": ("x", 6.0)},
    )

    # Surface selections may differ by rank. In particular, all ranks must participate
    # in boundary-ID discovery even when only rank 0 requests a boundary surface.
    rank_dependent_surface_prefix = "RankDependentSurfaces_p"
    rank_dependent_boundaries = ["xmin"] if rank == 0 else []
    rank_dependent_interior = {} if rank == 0 else {"inter_x": ("x", 6.0)}
    problem.WriteSurfaceAngularFluxes(
        rank_dependent_surface_prefix,
        boundary_surfaces=rank_dependent_boundaries,
        interior_surfaces=rank_dependent_interior,
    )

    # Every requested surface is readable on every rank, including ranks that own no
    # faces on that surface.
    boundary_surfaces = problem.ReadSurfaceAngularFluxes(
        forward_surface_prefix, ["xmin", "xmax", "ymin", "ymax"]
    )
    assert len(boundary_surfaces) == 4
    for surface in boundary_surfaces:
        if not surface["data"]["psi"]:
            if surface["data"]["node_index"] or surface["data"]["dir_index"]:
                raise RuntimeError("Empty surface contains angular-flux start offsets")

    # Adjoint solve and response evaluation.
    problem.SetAdjoint(True)
    adjoint_source = VolumetricSource(logical_volume=detector_volume, group_strength=[0.8])
    problem.SetVolumetricSources(volumetric_sources=[adjoint_source])
    solver.Execute()

    adjoint_flux_prefix = "AdjFluxMoments_p"
    adjoint_surface_prefix = "InteriorSurf_AdjSrc_p"
    problem.WriteFluxMoments(adjoint_flux_prefix)
    problem.WriteSurfaceAngularFluxes(
        adjoint_surface_prefix, interior_surfaces={"inter_x": ("x", 6.0)}
    )

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
    local_forward_flux = sum(forward_surface_up[0]["data"]["psi"], 0.0) + sum(
        forward_surface_down[0]["data"]["psi"], 0.0
    )
    forward_boundary_flux = global_sum(local_forward_flux)

    adjoint_surface_up = problem.ReadSurfaceAngularFluxes(
        adjoint_surface_prefix, ["inter_x_u"]
    )
    adjoint_surface_down = problem.ReadSurfaceAngularFluxes(
        adjoint_surface_prefix, ["inter_x_d"]
    )
    local_adjoint_importance = sum(adjoint_surface_up[0]["data"]["psi"], 0.0) + sum(
        adjoint_surface_down[0]["data"]["psi"], 0.0
    )
    adjoint_boundary_importance = global_sum(local_adjoint_importance)

    if rank == 0:
        print(f"Forward QoI={forward_qoi:.5e}")
        print(f"Adjoint QoI={adjoint_qoi:.5e}")
        print(f"Forward Boundary Flux={forward_boundary_flux:.5e}")
        print(f"Adjoint Boundary Importance={adjoint_boundary_importance:.5e}")

    barrier()
    for file_prefix in (
        forward_surface_prefix,
        rank_dependent_surface_prefix,
        adjoint_flux_prefix,
        adjoint_surface_prefix,
    ):
        os.remove(f"{file_prefix}{rank}.h5")
