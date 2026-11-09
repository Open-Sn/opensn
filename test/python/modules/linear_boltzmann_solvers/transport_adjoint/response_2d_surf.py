#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
2D transport response test exercising surface angular-flux I/O.

The interior plane x = 6 separates the source (x < 1) from the detector (x > 9.5) and all
exterior boundaries are vacuum, so the response is carried entirely across the plane:

    QoI = int_S int_4pi (Omega . n) psi(r, Omega) psi^dagger(r, Omega) dOmega dA,

with n = +x. For the upwind DFEM discretization the discrete identity is exact when psi is taken
from its upwind side of the face and psi^dagger from its own upwind side (the downstream side of
the forward direction). The surface QoI is evaluated from the exported directions, weights,
direction cosines, face mass matrices, and angular fluxes, and must equal the volumetric forward
QoI and the adjoint response. Evaluated along the inter_x_d normal, n = -x, which points from the
detector side to the source side, the surface response must equal minus the QoI.
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
    from pyopensn.mesh import KBAGraphPartitioner, OrthogonalMeshGenerator
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
    meshgen = OrthogonalMeshGenerator(
        node_sets=[nodes, nodes],
        partitioner=KBAGraphPartitioner(nx=1, ny=2, nz=1, ycuts=[length / 2.0]),
    )
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

    invalid_interior_cases = [
        ({"cuts_cells": ("x", 6.25)}, "cuts through mesh cells"),
        ({"outside": ("x", 12.0)}, "lies outside the mesh"),
        ({"first": ("y", 4.0), "second": ("y", 4.0)}, "defined on the same plane"),
    ]
    for interior_surfaces, expected_error in invalid_interior_cases:
        rejected = False
        try:
            problem.WriteSurfaceAngularFluxes(
                invalid_surface_prefix, interior_surfaces=interior_surfaces
            )
        except ValueError as error:
            rejected = expected_error in str(error)
        if not rejected:
            raise RuntimeError(f"Invalid interior surface {interior_surfaces} was not rejected")
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
            mapping = surface["mapping"]
            if mapping["normals"] or mapping["fe_shape"] or mapping["M_ij"]:
                raise RuntimeError("Empty surface contains face data")

    # Any file of a set can be read by index, so fewer ranks than files can reach every file.
    # Reading both partitions on each rank recovers every face of the x = 6 plane, one per cell row.
    (own_up,) = problem.ReadSurfaceAngularFluxes(forward_surface_prefix, ["inter_x_u"])
    (other_up,) = problem.ReadSurfaceAngularFluxes(
        forward_surface_prefix, ["inter_x_u"], partition=(rank + 1) % size
    )
    if own_up["num_partitions"] != size or other_up["num_partitions"] != size:
        raise RuntimeError("Unexpected number of partitions in the surface file set")
    if len(own_up["mapping"]["cell_ids"]) + len(other_up["mapping"]["cell_ids"]) != num_cells:
        raise RuntimeError("Reading every partition does not recover every face of the plane")

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
        data=[
            {"name": "detector", "file_prefixes": {"flux_moments": adjoint_flux_prefix}}
        ],
        sources={
            "material": [
                {"block_id": source_block_id, "strength": [source_strength / source_area]}
            ]
        },
    )
    adjoint_qoi = evaluator.EvaluateResponse("detector")
    try:
        evaluator.EvaluateResponse("missing")
        raise RuntimeError("An unknown data set was not rejected")
    except ValueError as error:
        if "No data named" not in str(error):
            raise

    def read_surface(prefix, tag, adjoint):
        (surface,) = problem.ReadSurfaceAngularFluxes(prefix, [tag])
        if surface["adjoint"] != adjoint:
            raise RuntimeError(f"Unexpected adjoint flag for {tag} in {prefix}")
        return surface

    def face_nodes(surface, face):
        mapping = surface["mapping"]
        start = sum(mapping["num_face_nodes"][:face])
        return [
            (mapping["nodes_x"][start + i], mapping["nodes_y"][start + i])
            for i in range(mapping["num_face_nodes"][face])
        ]

    def node_permutation(nodes_from, nodes_to):
        """Index in nodes_to of each node in nodes_from."""
        def distance_sq(k, x, y):
            return (nodes_to[k][0] - x) ** 2 + (nodes_to[k][1] - y) ** 2

        return [
            min(range(len(nodes_to)), key=lambda k: distance_sq(k, x, y)) for x, y in nodes_from
        ]

    forward_up = read_surface(forward_surface_prefix, "inter_x_u", False)
    forward_down = read_surface(forward_surface_prefix, "inter_x_d", False)
    adjoint_up = read_surface(adjoint_surface_prefix, "inter_x_u", True)
    adjoint_down = read_surface(adjoint_surface_prefix, "inter_x_d", True)

    # Both sides of each face must be local to evaluate the upwind traces. The partition used by
    # this test does not cut the x = 6 plane.
    if len(forward_up["mapping"]["cell_ids"]) != len(forward_down["mapping"]["cell_ids"]):
        raise RuntimeError("The x = 6 plane is split across ranks")

    # The groupset quadrature and the face normals give the direction cosines; one group.
    omega = forward_up["data"]["omega"]
    weights = forward_up["data"]["weights"]
    num_dirs = len(weights)
    local_surface_qoi = 0.0
    m_offset = 0
    for face, num_nodes in enumerate(forward_up["mapping"]["num_face_nodes"]):
        centroid_key = next(k for k, v in forward_up["mapping"]["cell_map"].items() if v == face)
        face_down = forward_down["mapping"]["cell_map"].get(centroid_key)
        if face_down is None:
            raise RuntimeError("The x = 6 plane is split across ranks")
        nodes_up = face_nodes(forward_up, face)
        down_of_up = node_permutation(nodes_up, face_nodes(forward_down, face_down))

        def psi(surface, face_index, node, d):
            k = sum(surface["mapping"]["num_face_nodes"][:face_index]) + node
            return surface["data"]["psi"][k * num_dirs + d]

        normal = forward_up["mapping"]["normals"][3 * face : 3 * face + 3]
        mass = forward_up["mapping"]["M_ij"][m_offset : m_offset + num_nodes * num_nodes]
        for d in range(num_dirs):
            mu = sum(omega[3 * d + c] * normal[c] for c in range(3))
            weight = weights[d]
            for i in range(num_nodes):
                if mu > 0.0:
                    psi_i = psi(forward_up, face, i, d)
                else:
                    psi_i = psi(forward_down, face_down, down_of_up[i], d)
                for j in range(num_nodes):
                    if mu > 0.0:
                        psi_adj_j = psi(adjoint_down, face_down, down_of_up[j], d)
                    else:
                        psi_adj_j = psi(adjoint_up, face, j, d)
                    local_surface_qoi += weight * mu * mass[i * num_nodes + j] * psi_i * psi_adj_j
        m_offset += num_nodes * num_nodes
    surface_qoi = global_sum(local_surface_qoi)

    # The evaluator matches faces by location and redistributes them over the ranks, so it does not
    # require both sides of the plane on one rank.
    # The interior surface name reads both stored sides, inter_x_u and inter_x_d.
    interior_tags = ["inter_x"]
    evaluator.SetOptions(
        data=[
            {
                "name": "forward_all",
                "file_prefixes": {"surface_angular_fluxes": forward_surface_prefix},
                "surfaces": ["xmin", "xmax", "ymin", "ymax"] + interior_tags,
            },
            {
                "name": "adjoint_interior",
                "file_prefixes": {"surface_angular_fluxes": adjoint_surface_prefix},
                "surfaces": interior_tags,
            },
        ]
    )
    evaluator_surface = evaluator.EvaluateSurfaceResponse(
        "forward_all", "adjoint_interior", surfaces=["inter_x"]
    )["inter_x"]
    evaluator_surface_qoi = evaluator_surface["total"]
    # The side tag inter_x_d measures the same interior surface along its own normal, n = -x.
    reversed_surface = evaluator.EvaluateSurfaceResponse(
        "forward_all", "adjoint_interior", surfaces=["inter_x_d"]
    )["inter_x_d"]
    reversed_surface_qoi = reversed_surface["total"]
    # Along -x the parts swap and change sign, and each total is the sum of its parts.
    if (
        reversed_surface["outgoing"] != -evaluator_surface["incoming"]
        or reversed_surface["incoming"] != -evaluator_surface["outgoing"]
        or evaluator_surface["total"]
        != evaluator_surface["outgoing"] + evaluator_surface["incoming"]
    ):
        raise RuntimeError("Inconsistent outgoing and incoming parts of inter_x")

    # By default every forward surface is evaluated; the adjoint data lack the boundary tags.
    try:
        evaluator.EvaluateSurfaceResponse("forward_all", "adjoint_interior")
        raise RuntimeError("A surface missing from the adjoint data was not rejected")
    except ValueError as error:
        if "requires tag" not in str(error):
            raise
    try:
        evaluator.EvaluateSurfaceResponse("adjoint_interior", "forward_all", surfaces=["inter_x"])
        raise RuntimeError("Swapped forward and adjoint data were not rejected")
    except ValueError as error:
        if "adjoint mode" not in str(error):
            raise

    if rank == 0:
        print(f"Forward QoI={forward_qoi:.8e}")
        print(f"Adjoint QoI={adjoint_qoi:.8e}")
        print(f"Surface QoI={surface_qoi:.8e}")
        print(f"Evaluator QoI={evaluator_surface_qoi:.8e}")
        # Relative errors of the surface identity: each surface QoI must equal the forward QoI.
        print(f"Reference identity error={abs(surface_qoi / forward_qoi - 1.0):.3e}")
        print(f"Evaluator identity error={abs(evaluator_surface_qoi / forward_qoi - 1.0):.3e}")
        # Along n = -x the surface QoI must equal minus the forward QoI.
        print(f"Reversed-normal identity error={abs(reversed_surface_qoi / forward_qoi + 1.0):.3e}")

    barrier()
    for file_prefix in (
        forward_surface_prefix,
        rank_dependent_surface_prefix,
        adjoint_flux_prefix,
        adjoint_surface_prefix,
    ):
        os.remove(f"{file_prefix}{rank}.h5")
