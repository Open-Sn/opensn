#!/usr/bin/env python3
"""Verify 3D multigroup uncollided attenuation.

Uncollided groups are independent, so the two-group solution must equal two
one-group solutions with each group's total cross section and strength to
round-off. Each group is also compared with the exact point-source solution
within the spatial discretization error of the mesh.
"""

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
uncollided_3d = uncollided_utils.uncollided_3d


def solve(grid, xs_map, num_groups, source, strength, sample_points, file_name):
    # Small region around the source
    near_source_region = RPPLogicalVolume(
        xmin=source[0] - 0.15,
        xmax=source[0] + 0.15,
        ymin=source[1] - 0.15,
        ymax=source[1] + 0.15,
        zmin=source[2] - 0.15,
        zmax=source[2] + 0.15,
    )

    remove_file(file_name)
    problem = UncollidedProblem(
        mesh=grid,
        num_groups=num_groups,
        groupsets=[{"groups_from_to": [0, num_groups - 1]}],
        xs_map=xs_map,
        point_sources=[PointSource(location=list(source), strength=strength)],
        near_source=[near_source_region],
        scattering_order=0,
    )
    solver = UncollidedSolver(problem=problem, file_name=file_name)
    solver.Initialize()
    solver.Execute()

    scalar_fluxes = problem.GetScalarFluxFieldFunction()
    values = [
        [
            point_value(scalar_fluxes[group], point, FieldFunctionInterpolationPoint, Vector3)
            for point in sample_points
        ]
        for group in range(num_groups)
    ]
    remove_file(file_name)
    return values


if __name__ == "__main__":
    if size != 1:
        sys.exit(f"Expected one process, got {size}.")

    grid = extruded_mesh(
        ExtruderMeshGenerator, FromFileMeshGenerator, "triangle_mesh2x2.obj", num_layers=16
    )
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

    source = (0.11, -0.23, 1.07)
    strength = [1.0, 0.35]
    sample_points = [
        (-0.42, 0.16, 0.81),
        (0.36, 0.29, 1.33),
        (0.51, -0.34, 0.92),
    ]

    multigroup = solve(
        grid,
        [{"block_ids": [0], "xs": xs}],
        2,
        source,
        strength,
        sample_points,
        "uncollided_3d_multigroup_analytic.h5",
    )

    max_group_error = 0.0
    max_decoupling_difference = 0.0
    for group in range(2):
        group_xs = MultiGroupXS()
        group_xs.CreateSimpleOneGroup(sigma_t=sigma_t[group], c=0.0)
        single = solve(
            grid,
            [{"block_ids": [0], "xs": group_xs}],
            1,
            source,
            [strength[group]],
            sample_points,
            f"uncollided_3d_multigroup_group{group}.h5",
        )[0]
        for point, numeric, single_value in zip(sample_points, multigroup[group], single):
            exact = uncollided_3d(strength[group], sigma_t[group], source, point)
            max_group_error = max(max_group_error, relative_error(numeric, exact))
            max_decoupling_difference = max(
                max_decoupling_difference, relative_error(numeric, single_value)
            )

    if rank == 0:
        print(f"Uncollided3DMultigroupDecouplingDifference={max_decoupling_difference:.8e}")
        print(f"Uncollided3DMultigroupMaxError={max_group_error:.8e}")
        sys.stdout.flush()

    if max_decoupling_difference > 1.0e-12:
        raise RuntimeError(
            f"Multigroup flux differs from one-group solutions: {max_decoupling_difference}"
        )
    if max_group_error > 0.08:
        raise RuntimeError(f"3D multigroup scalar-flux error is too large: {max_group_error}")
