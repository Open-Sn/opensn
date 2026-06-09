#!/usr/bin/env python3
"""Optical-path attenuation through oblique material interfaces."""

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
    from pyopensn.mesh import FromFileMeshGenerator
    from pyopensn.solver import UncollidedProblem, UncollidedSolver
    from pyopensn.source import PointSource
    from pyopensn.xs import MultiGroupXS

from uncollided_unstructured_utils import mesh_path, point_value, relative_error, remove_file


def reference_flux(source, point, sigma_t):
    source_u = source[0] + 0.35 * source[1]
    point_u = point[0] + 0.35 * point[1]
    distance = math.dist(source[:2], point[:2])
    first_fraction = (-0.25 - source_u) / (point_u - source_u)
    second_fraction = (0.25 - source_u) / (point_u - source_u)
    lengths = (
        distance * first_fraction,
        distance * (second_fraction - first_fraction),
        distance * (1.0 - second_fraction),
    )
    optical_path = sum(xs * length for xs, length in zip(sigma_t, lengths))
    return math.exp(-optical_path) / (2.0 * math.pi * distance)


if __name__ == "__main__":
    if size != 1:
        sys.exit(f"Expected one process, got {size}.")

    grid = FromFileMeshGenerator(
        filename=mesh_path("uncollided_oblique_layers.msh")
    ).Execute()
    sigma_t = (0.2, 0.8, 1.4)
    cross_sections = []
    for value in sigma_t:
        xs = MultiGroupXS()
        xs.CreateSimpleOneGroup(sigma_t=value, c=0.0)
        cross_sections.append(xs)

    source = (-0.8, -0.4, 0.0)
    whole_domain = RPPLogicalVolume(
        xmin=-1.01,
        xmax=1.01,
        ymin=-1.01,
        ymax=1.01,
        infz=True,
    )
    file_name = "uncollided_2d_oblique_layers.h5"
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
        near_source=[whole_domain],
        scattering_order=0,
    )
    solver = UncollidedSolver(problem=problem, file_name=file_name)
    solver.Initialize()
    solver.Execute()

    scalar_flux = problem.GetScalarFluxFieldFunction()[0]
    sample_points = [
        (0.75, 0.40, 0.0),
        (0.60, 0.70, 0.0),
        (0.80, 0.00, 0.0),
    ]
    max_error = max(
        relative_error(
            point_value(scalar_flux, point, FieldFunctionInterpolationPoint, Vector3),
            reference_flux(source, point, sigma_t),
        )
        for point in sample_points
    )

    if rank == 0:
        print(f"Uncollided2DObliqueLayersError={max_error:.8e}")

    remove_file(file_name)
    if max_error > 0.02:
        raise RuntimeError(f"Oblique-interface optical-path error is too large: {max_error}")
