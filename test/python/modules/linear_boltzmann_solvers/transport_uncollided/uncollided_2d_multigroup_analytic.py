#!/usr/bin/env python3
"""Verify 2D multigroup uncollided attenuation."""

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
relative_error = uncollided_utils.relative_error
remove_file = uncollided_utils.remove_file
uncollided_2d = uncollided_utils.uncollided_2d


if __name__ == "__main__":
    if size != 1:
        sys.exit(f"Expected one process, got {size}.")

    grid = FromFileMeshGenerator(filename=mesh_path("triangle_mesh2x2_fine.obj")).Execute()
    grid.SetUniformBlockID(0)

    xs = MultiGroupXS()
    xs.LoadFromOpenSn(
        os.path.abspath(
            os.path.join(
                os.path.dirname(__file__),
                "../../../../assets/xs/simple_2g_downscatter_td.cxs",
            )
        )
    )
    sigma_t = xs.sigma_t

    source = (0.11, -0.23, 0.0)
    strength = [1.0, 0.35]
    sample_points = [
        (-0.42, 0.16, 0.0),
        (0.36, 0.29, 0.0),
        (0.51, -0.34, 0.0),
    ]

    whole_domain = RPPLogicalVolume(
        xmin=-1.01,
        xmax=1.01,
        ymin=-1.01,
        ymax=1.01,
        infz=True,
    )

    file_name = "uncollided_2d_multigroup_analytic.h5"
    remove_file(file_name)
    problem = UncollidedProblem(
        mesh=grid,
        num_groups=2,
        groupsets=[{"groups_from_to": [0, 1]}],
        xs_map=[{"block_ids": [0], "xs": xs}],
        point_sources=[PointSource(location=list(source), strength=strength)],
        near_source=[whole_domain],
        scattering_order=0,
    )
    solver = UncollidedSolver(problem=problem, file_name=file_name)
    solver.Initialize()
    solver.Execute()

    scalar_fluxes = problem.GetScalarFluxFieldFunction()
    max_group_error = 0.0
    for group in range(2):
        for point in sample_points:
            numeric = point_value(
                scalar_fluxes[group], point, FieldFunctionInterpolationPoint, Vector3
            )
            exact = uncollided_2d(strength[group], sigma_t[group], source, point)
            max_group_error = max(max_group_error, relative_error(numeric, exact))

    remove_file(file_name)

    if rank == 0:
        print(f"Uncollided2DMultigroupMaxError={max_group_error:.8e}")
        sys.stdout.flush()

    if max_group_error > 1.0e-3:
        raise RuntimeError(f"2D multigroup scalar-flux error is too large: {max_group_error}")
