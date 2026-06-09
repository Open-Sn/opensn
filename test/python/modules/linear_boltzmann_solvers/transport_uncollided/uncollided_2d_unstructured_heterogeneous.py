#!/usr/bin/env python3
"""Heterogeneous 2D optical-path attenuation on an unstructured mesh."""

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


if __name__ == "__main__":
    if size != 1:
        sys.exit(f"Expected one process, got {size}.")

    grid = FromFileMeshGenerator(filename=mesh_path("triangle_mesh2x2_fine.obj")).Execute()
    grid.SetUniformBlockID(0)
    right_half = RPPLogicalVolume(xmin=0.0, xmax=1.01, infy=True, infz=True)
    grid.SetBlockIDFromLogicalVolume(right_half, 1, True)

    sigma_left = 0.25
    sigma_right = 1.15
    xs_left = MultiGroupXS()
    xs_left.CreateSimpleOneGroup(sigma_t=sigma_left, c=0.0)
    xs_right = MultiGroupXS()
    xs_right.CreateSimpleOneGroup(sigma_t=sigma_right, c=0.0)

    source = (-0.62, 0.13, 0.0)
    point = (0.61, 0.13, 0.0)
    point_source = PointSource(location=list(source), strength=[1.0])
    whole_domain = RPPLogicalVolume(
        xmin=-1.01,
        xmax=1.01,
        ymin=-1.01,
        ymax=1.01,
        infz=True,
    )

    file_name = "uncollided_2d_unstructured_heterogeneous.h5"
    remove_file(file_name)
    problem = UncollidedProblem(
        mesh=grid,
        num_groups=1,
        groupsets=[{"groups_from_to": [0, 0]}],
        xs_map=[
            {"block_ids": [0], "xs": xs_left},
            {"block_ids": [1], "xs": xs_right},
        ],
        point_sources=[point_source],
        near_source=[whole_domain],
        scattering_order=0,
    )
    solver = UncollidedSolver(problem=problem, file_name=file_name)
    solver.Initialize()
    solver.Execute()

    scalar_flux = problem.GetScalarFluxFieldFunction()[0]
    value = point_value(scalar_flux, point, FieldFunctionInterpolationPoint, Vector3)
    radius = point[0] - source[0]
    optical_path = sigma_left * (0.0 - source[0]) + sigma_right * point[0]
    reference = math.exp(-optical_path) / (2.0 * math.pi * radius)
    error = relative_error(value, reference)

    if rank == 0:
        print(f"Uncollided2DHeterogeneousValue={value:.8e}")
        print(f"Uncollided2DHeterogeneousReference={reference:.8e}")
        print(f"Uncollided2DHeterogeneousRelativeError={error:.8e}")

    remove_file(file_name)
    if error > 0.02:
        raise RuntimeError(f"Heterogeneous optical-path error is too large: {error}")
