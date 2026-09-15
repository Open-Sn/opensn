#!/usr/bin/env python3
"""Homogeneous 3D tetrahedral-mesh convergence test."""

import importlib
import os
import sys

if "opensn_console" not in globals():
    from mpi4py import MPI

    size = MPI.COMM_WORLD.size
    rank = MPI.COMM_WORLD.rank
    sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), "../../../../../")))
    from pyopensn.fieldfunc import (
        FieldFunctionInterpolationPoint,
    )
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
uncollided_3d = uncollided_utils.uncollided_3d


def compute_error(mesh_name, file_name):
    grid = FromFileMeshGenerator(filename=mesh_path(mesh_name)).Execute()
    grid.SetUniformBlockID(0)

    sigma_t = 35.0
    xs = MultiGroupXS()
    xs.CreateSimpleOneGroup(sigma_t=sigma_t, c=0.0)

    # Chosen to clear the tetrahedral meshes' cell faces on all three
    # resolutions. 3D meshes leave less margin than 2D, so this may still
    # trigger a warning.
    source = (0.01758, 0.01118, 0.012957)
    point_source = PointSource(location=list(source), strength=[1.0])
    # Small region around the source
    near_source_region = RPPLogicalVolume(
        xmin=source[0] - 0.008,
        xmax=source[0] + 0.008,
        ymin=source[1] - 0.008,
        ymax=source[1] + 0.008,
        zmin=source[2] - 0.008,
        zmax=source[2] + 0.008,
    )

    remove_file(file_name)
    problem = UncollidedProblem(
        mesh=grid,
        num_groups=1,
        groupsets=[{"groups_from_to": [0, 0]}],
        xs_map=[{"block_ids": [0], "xs": xs}],
        point_sources=[point_source],
        near_source=[near_source_region],
        scattering_order=0,
    )
    solver = UncollidedSolver(problem=problem, file_name=file_name)
    solver.Initialize()
    solver.Execute()

    scalar_flux = problem.GetScalarFluxFieldFunction()[0]
    sample_points = [
        (0.0230, 0.0153, 0.0161),
        (0.0147, 0.0240, 0.0170),
        (0.0090, 0.0110, 0.0230),
    ]
    point_errors = []
    for point in sample_points:
        value = point_value(scalar_flux, point, FieldFunctionInterpolationPoint, Vector3)
        reference = uncollided_3d(1.0, sigma_t, source, point)
        point_errors.append(relative_error(value, reference))

    remove_file(file_name)
    return point_errors


if __name__ == "__main__":
    if size != 1:
        sys.exit(f"Expected one process, got {size}.")

    meshes = [
        "uncollided_cube_coarse.msh",
        "cube3.2.msh",
        "uncollided_cube_fine.msh",
    ]
    results = [
        compute_error(mesh_name, f"uncollided_3d_convergence_{level}.h5")
        for level, mesh_name in enumerate(meshes)
    ]
    point_errors = results
    errors = [
        sum(error * error for error in mesh_errors) ** 0.5 / len(mesh_errors) ** 0.5
        for mesh_errors in point_errors
    ]
    max_errors = [max(mesh_errors) for mesh_errors in point_errors]

    if rank == 0:
        print(f"Uncollided3DCoarseError={errors[0]:.8e}")
        print(f"Uncollided3DMediumError={errors[1]:.8e}")
        print(f"Uncollided3DFineError={errors[2]:.8e}")
        print(f"Uncollided3DFineMaxError={max_errors[2]:.8e}")

    if not errors[1] < 0.8 * errors[0]:
        raise RuntimeError(f"3D coarse-to-medium convergence failed: {errors}")
    if not errors[2] < 0.9 * errors[1]:
        raise RuntimeError(f"3D medium-to-fine convergence failed: {errors}")
    if max_errors[2] > 0.11:
        raise RuntimeError(f"3D fine-mesh scalar-flux error is too large: {max_errors[2]}")
