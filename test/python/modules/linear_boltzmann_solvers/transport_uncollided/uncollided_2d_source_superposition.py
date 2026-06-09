#!/usr/bin/env python3
"""Verify linear superposition for multiple 2D point sources."""

import importlib
import os
import sys

if "opensn_console" not in globals():
    from mpi4py import MPI

    rank = MPI.COMM_WORLD.rank
    size = MPI.COMM_WORLD.size
    sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), "../../../../../")))
    from pyopensn.fieldfunc import FieldFunctionInterpolationPoint
    from pyopensn.logvol import RPPLogicalVolume
    from pyopensn.math import Vector3
    from pyopensn.mesh import FromFileMeshGenerator
    from pyopensn.solver import UncollidedProblem, UncollidedSolver
    from pyopensn.source import PointSource
    from pyopensn.xs import MultiGroupXS

sys.path.append(os.path.dirname(__file__))
uncollided_utils = importlib.import_module("uncollided_unstructured_utils")
mesh_path = uncollided_utils.mesh_path
point_value = uncollided_utils.point_value
remove_file = uncollided_utils.remove_file


def solve_scalar_flux(grid, xs, whole_domain, sources, file_name):
    remove_file(file_name)
    problem = UncollidedProblem(
        mesh=grid,
        num_groups=1,
        groupsets=[{"groups_from_to": [0, 0]}],
        xs_map=[{"block_ids": [0], "xs": xs}],
        point_sources=[
            PointSource(location=list(source), strength=[strength])
            for source, strength in sources
        ],
        near_source=[whole_domain] * len(sources),
        scattering_order=0,
    )
    solver = UncollidedSolver(problem=problem, file_name=file_name)
    solver.Initialize()
    solver.Execute()
    scalar_flux = problem.GetScalarFluxFieldFunction()[0]
    return scalar_flux


if __name__ == "__main__":
    if size != 1:
        sys.exit(f"Expected one process, got {size}.")

    grid = FromFileMeshGenerator(filename=mesh_path("triangle_mesh2x2_fine.obj")).Execute()
    grid.SetUniformBlockID(0)

    xs = MultiGroupXS()
    xs.CreateSimpleOneGroup(sigma_t=0.6, c=0.0)

    whole_domain = RPPLogicalVolume(
        xmin=-1.01,
        xmax=1.01,
        ymin=-1.01,
        ymax=1.01,
        infz=True,
    )

    source_a = ((-0.31, 0.18, 0.0), 1.0)
    source_b = ((0.29, -0.27, 0.0), 0.7)
    sample_points = [
        (-0.52, -0.41, 0.0),
        (-0.18, 0.49, 0.0),
        (0.47, -0.08, 0.0),
        (0.14, 0.34, 0.0),
    ]

    scalar_a = solve_scalar_flux(
        grid, xs, whole_domain, [source_a], "uncollided_superposition_a.h5"
    )
    scalar_b = solve_scalar_flux(
        grid, xs, whole_domain, [source_b], "uncollided_superposition_b.h5"
    )
    scalar_ab = solve_scalar_flux(
        grid,
        xs,
        whole_domain,
        [source_a, source_b],
        "uncollided_superposition_ab.h5",
    )

    max_relative_difference = 0.0
    for point in sample_points:
        value_a = point_value(scalar_a, point, FieldFunctionInterpolationPoint, Vector3)
        value_b = point_value(scalar_b, point, FieldFunctionInterpolationPoint, Vector3)
        value_ab = point_value(scalar_ab, point, FieldFunctionInterpolationPoint, Vector3)
        reference = max(abs(value_a + value_b), 1.0e-14)
        max_relative_difference = max(
            max_relative_difference,
            abs(value_ab - (value_a + value_b)) / reference,
        )

    remove_file("uncollided_superposition_a.h5")
    remove_file("uncollided_superposition_b.h5")
    remove_file("uncollided_superposition_ab.h5")

    if rank == 0:
        print(f"Uncollided2DSuperpositionRelativeDifference={max_relative_difference:.8e}")
        sys.stdout.flush()

    if max_relative_difference > 1.0e-12:
        raise RuntimeError(
            f"Multiple point-source solution broke linear superposition: {max_relative_difference}"
        )
