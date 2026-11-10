#!/usr/bin/env python3
"""Optical-path attenuation through oblique material interfaces in 3D.

The 2D oblique-layer mesh is extruded in z, so the interfaces are planes
x + 0.35 y = const. The layered problem is compared with a homogeneous problem
on the same mesh whose cross section equals that of the source layer. The exact
ratio of the two fluxes at a point is exp(-(tau_layered - tau_homogeneous)), and
most of the spatial discretization error cancels in the computed ratio, so the
ratio checks the interface attenuation tightly. The absolute error is limited by
the mesh resolution and is checked more loosely.
"""

import importlib
import math
import os
import sys

if "opensn_console" not in globals():
    from mpi4py import MPI

    size = MPI.COMM_WORLD.size
    rank = MPI.COMM_WORLD.rank
    sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), "../../../../../")))
    from pyopensn.fieldfunc import FieldFunctionInterpolationPoint
    from pyopensn.logvol import RPPLogicalVolume
    from pyopensn.math import Vector3
    from pyopensn.mesh import ExtruderMeshGenerator, FromFileMeshGenerator
    from pyopensn.solver import UncollidedProblem, UncollidedSolver
    from pyopensn.source import PointSource
    from pyopensn.xs import MultiGroupXS

sys.path.append(os.path.dirname(__file__))
uncollided_utils = importlib.import_module("uncollided_unstructured_utils")
extruded_mesh = uncollided_utils.extruded_mesh
point_value = uncollided_utils.point_value
relative_error = uncollided_utils.relative_error
remove_file = uncollided_utils.remove_file


def reference_flux(source, point, sigma_t):
    source_u = source[0] + 0.35 * source[1]
    point_u = point[0] + 0.35 * point[1]
    distance = math.dist(source, point)
    first_fraction = (-0.25 - source_u) / (point_u - source_u)
    second_fraction = (0.25 - source_u) / (point_u - source_u)
    lengths = (
        distance * first_fraction,
        distance * (second_fraction - first_fraction),
        distance * (1.0 - second_fraction),
    )
    optical_path = sum(xs * length for xs, length in zip(sigma_t, lengths))
    return math.exp(-optical_path) / (4.0 * math.pi * distance**2)


def solve(grid, sigma_t, source, sample_points, file_name):
    cross_sections = []
    for value in sigma_t:
        xs = MultiGroupXS()
        xs.CreateSimpleOneGroup(sigma_t=value, c=0.0)
        cross_sections.append(xs)

    # Small region around the source
    near_source_region = RPPLogicalVolume(
        xmin=source[0] - 0.08,
        xmax=source[0] + 0.08,
        ymin=source[1] - 0.08,
        ymax=source[1] + 0.08,
        zmin=source[2] - 0.08,
        zmax=source[2] + 0.08,
    )
    remove_file(file_name)
    problem = UncollidedProblem(
        mesh=grid,
        num_groups=1,
        groupsets=[{"groups_from_to": [0, 0]}],
        xs_map=[
            {"block_ids": [block_id], "xs": xs}
            for block_id, xs in enumerate(cross_sections, start=1)
        ],
        point_sources=[PointSource(location=list(source), strength=[1.0])],
        near_source=[near_source_region],
        scattering_order=0,
    )
    solver = UncollidedSolver(problem=problem, file_name=file_name)
    solver.Initialize()
    solver.Execute()

    scalar_flux = problem.GetScalarFluxFieldFunction()[0]
    values = [
        point_value(scalar_flux, point, FieldFunctionInterpolationPoint, Vector3)
        for point in sample_points
    ]
    remove_file(file_name)
    return values


if __name__ == "__main__":
    if size != 1:
        sys.exit(f"Expected one process, got {size}.")

    grid = extruded_mesh(
        ExtruderMeshGenerator,
        FromFileMeshGenerator,
        "uncollided_oblique_layers.msh",
        height=1.0,
        num_layers=8,
    )
    layered = (0.2, 0.8, 1.4)
    homogeneous = (0.2, 0.2, 0.2)

    # Within the first (u < -0.25) material layer
    source = (-0.778517, -0.390592, 0.53)
    sample_points = [
        (0.75, 0.40, 0.71),
        (0.60, 0.70, 0.36),
        (0.80, 0.00, 0.62),
    ]
    layered_values = solve(
        grid, layered, source, sample_points, "uncollided_3d_oblique_layers.h5"
    )
    homogeneous_values = solve(
        grid, homogeneous, source, sample_points, "uncollided_3d_oblique_homogeneous.h5"
    )

    max_error = 0.0
    max_ratio_error = 0.0
    for point, layered_value, homogeneous_value in zip(
        sample_points, layered_values, homogeneous_values
    ):
        layered_reference = reference_flux(source, point, layered)
        homogeneous_reference = reference_flux(source, point, homogeneous)
        max_error = max(max_error, relative_error(layered_value, layered_reference))
        max_ratio_error = max(
            max_ratio_error,
            relative_error(
                layered_value / homogeneous_value, layered_reference / homogeneous_reference
            ),
        )

    if rank == 0:
        print(f"Uncollided3DObliqueLayersRatioError={max_ratio_error:.8e}")
        print(f"Uncollided3DObliqueLayersError={max_error:.8e}")

    if max_ratio_error > 2.0e-3:
        raise RuntimeError(f"Oblique-interface attenuation error is too large: {max_ratio_error}")
    if max_error > 0.15:
        raise RuntimeError(f"Oblique-interface optical-path error is too large: {max_error}")
