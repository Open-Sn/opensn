#!/usr/bin/env python3
"""Heterogeneous 3D optical-path attenuation across a planar material interface.

The orthogonal mesh has a node plane at x = 0, so the interface follows cell faces.
The two-material problem is compared with a homogeneous problem on the same mesh
whose cross section equals that of the source side. The exact ratio of the two
fluxes at a point is exp(-(tau_heterogeneous - tau_homogeneous)), and most of the
spatial discretization error cancels in the computed ratio, so the ratio checks
the interface attenuation tightly. The absolute error is limited by the mesh
resolution and is checked more loosely.
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
    from pyopensn.mesh import OrthogonalMeshGenerator
    from pyopensn.solver import UncollidedProblem, UncollidedSolver
    from pyopensn.source import PointSource
    from pyopensn.xs import MultiGroupXS

sys.path.append(os.path.dirname(__file__))
uncollided_utils = importlib.import_module("uncollided_unstructured_utils")
point_value = uncollided_utils.point_value
relative_error = uncollided_utils.relative_error
remove_file = uncollided_utils.remove_file


def optical_path(source, point, sigma_left, sigma_right):
    # The segment crosses the x = 0 interface; split its length by the x fractions.
    left_fraction = (0.0 - source[0]) / (point[0] - source[0])
    return math.dist(source, point) * (
        sigma_left * left_fraction + sigma_right * (1.0 - left_fraction)
    )


def solve(grid, sigma_left, sigma_right, source, point, file_name):
    xs_left = MultiGroupXS()
    xs_left.CreateSimpleOneGroup(sigma_t=sigma_left, c=0.0)
    xs_right = MultiGroupXS()
    xs_right.CreateSimpleOneGroup(sigma_t=sigma_right, c=0.0)

    # Small region around the source kept within the left material so it doesn't
    # cross the x=0 interface
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
            {"block_ids": [0], "xs": xs_left},
            {"block_ids": [1], "xs": xs_right},
        ],
        point_sources=[PointSource(location=list(source), strength=[1.0])],
        near_source=[near_source_region],
        scattering_order=0,
    )
    solver = UncollidedSolver(problem=problem, file_name=file_name)
    solver.Initialize()
    solver.Execute()

    scalar_flux = problem.GetScalarFluxFieldFunction()[0]
    value = point_value(scalar_flux, point, FieldFunctionInterpolationPoint, Vector3)
    remove_file(file_name)
    return value


if __name__ == "__main__":
    if size != 1:
        sys.exit(f"Expected one process, got {size}.")

    nodes = [-1.0 + 0.1 * i for i in range(21)]
    grid = OrthogonalMeshGenerator(node_sets=[nodes, nodes, nodes]).Execute()
    grid.SetUniformBlockID(0)
    right_half = RPPLogicalVolume(xmin=0.0, xmax=1.01, infy=True, infz=True)
    grid.SetBlockIDFromLogicalVolume(right_half, 1, True)

    sigma_left = 0.25
    sigma_right = 1.15
    # Left of the x=0 material interface
    source = (-0.589726, 0.133007, 0.031)
    point = (0.61, 0.13, -0.27)

    value = solve(grid, sigma_left, sigma_right, source, point, "uncollided_3d_het.h5")
    homogeneous_value = solve(grid, sigma_left, sigma_left, source, point, "uncollided_3d_hom.h5")

    radius = math.dist(source, point)
    tau = optical_path(source, point, sigma_left, sigma_right)
    reference = math.exp(-tau) / (4.0 * math.pi * radius**2)
    error = relative_error(value, reference)
    reference_ratio = math.exp(-(tau - optical_path(source, point, sigma_left, sigma_left)))
    ratio_error = relative_error(value / homogeneous_value, reference_ratio)

    if rank == 0:
        print(f"Uncollided3DHeterogeneousRatioError={ratio_error:.8e}")
        print(f"Uncollided3DHeterogeneousRelativeError={error:.8e}")

    if ratio_error > 1.0e-2:
        raise RuntimeError(f"Heterogeneous attenuation error is too large: {ratio_error}")
    if error > 0.1:
        raise RuntimeError(f"Heterogeneous optical-path error is too large: {error}")
