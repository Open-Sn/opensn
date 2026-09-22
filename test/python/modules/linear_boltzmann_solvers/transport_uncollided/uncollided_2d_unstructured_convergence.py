#!/usr/bin/env python3
"""Homogeneous 2D accuracy test on several independent unstructured meshes."""

import importlib
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

sys.path.append(os.path.dirname(__file__))
uncollided_utils = importlib.import_module("uncollided_unstructured_utils")
mesh_path = uncollided_utils.mesh_path
point_value = uncollided_utils.point_value
relative_error = uncollided_utils.relative_error
remove_file = uncollided_utils.remove_file
uncollided_2d = uncollided_utils.uncollided_2d


def compute_error(mesh_name, file_name):
    grid = FromFileMeshGenerator(filename=mesh_path(mesh_name)).Execute()
    grid.SetUniformBlockID(0)

    sigma_t = 0.7
    xs = MultiGroupXS()
    xs.CreateSimpleOneGroup(sigma_t=sigma_t, c=0.0)
    # Chosen so it isn't flush against a cell face on any of the three
    # meshes
    source = (0.04735, -0.033933, 0.0)
    near_source_region = RPPLogicalVolume(
        xmin=source[0] - 0.08,
        xmax=source[0] + 0.08,
        ymin=source[1] - 0.08,
        ymax=source[1] + 0.08,
        infz=True,
    )

    remove_file(file_name)
    problem = UncollidedProblem(
        mesh=grid,
        num_groups=1,
        groupsets=[{"groups_from_to": [0, 0]}],
        xs_map=[{"block_ids": [0], "xs": xs}],
        point_sources=[PointSource(location=list(source), strength=[1.0])],
        near_source=[near_source_region],
        scattering_order=0,
    )
    solver = UncollidedSolver(problem=problem, file_name=file_name)
    solver.Initialize()
    solver.Execute()

    scalar_flux = problem.GetScalarFluxFieldFunction()[0]
    sample_points = [
        (0.46, 0.11, 0.0),
        (-0.39, 0.43, 0.0),
        (0.22, -0.58, 0.0),
    ]
    error = max(
        relative_error(
            point_value(scalar_flux, point, FieldFunctionInterpolationPoint, Vector3),
            uncollided_2d(1.0, sigma_t, source, point),
        )
        for point in sample_points
    )
    remove_file(file_name)
    return error


if __name__ == "__main__":
    if size != 1:
        sys.exit(f"Expected one process, got {size}.")

    meshes = [
        "triangle_mesh2x2.obj",
        "triangle_mesh2x2_fine.obj",
        "triangle_mesh2x2_super_fine.obj",
    ]
    errors = [
        compute_error(mesh_name, f"uncollided_2d_convergence_{level}.h5")
        for level, mesh_name in enumerate(meshes)
    ]

    if rank == 0:
        print(f"Uncollided2DCoarseError={errors[0]:.8e}")
        print(f"Uncollided2DMediumError={errors[1]:.8e}")
        print(f"Uncollided2DFineError={errors[2]:.8e}")

    # These independently generated meshes do not form a nested refinement
    # sequence, so each result is checked against the analytic solution.
    if errors[0] > 0.06 or errors[1] > 0.20 or errors[2] > 0.15:
        raise RuntimeError(f"2D unstructured-mesh error is too large: {errors}")
