#!/usr/bin/env python3
"""Homogeneous 3D tetrahedral-mesh convergence test."""

import os
import sys

if "opensn_console" not in globals():
    from mpi4py import MPI

    size = MPI.COMM_WORLD.size
    rank = MPI.COMM_WORLD.rank
    sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), "../../../../../")))
    from pyopensn.fieldfunc import (
        FieldFunctionInterpolationPoint,
        FieldFunctionInterpolationVolume,
    )
    from pyopensn.logvol import RPPLogicalVolume
    from pyopensn.math import Vector3
    from pyopensn.mesh import FromFileMeshGenerator
    from pyopensn.solver import UncollidedProblem
    from pyopensn.source import PointSource
    from pyopensn.xs import MultiGroupXS

from uncollided_unstructured_utils import (
    mesh_path,
    point_value,
    relative_error,
    remove_file,
    uncollided_3d,
    volume_minimum,
)


def compute_error(mesh_name, file_name):
    grid = FromFileMeshGenerator(filename=mesh_path(mesh_name)).Execute()
    grid.SetUniformBlockID(0)

    sigma_t = 35.0
    xs = MultiGroupXS()
    xs.CreateSimpleOneGroup(sigma_t=sigma_t, c=0.0)

    source = (0.0147, 0.0153, 0.0161)
    point_source = PointSource(location=list(source), strength=[1.0])
    whole_domain = RPPLogicalVolume(
        xmin=-1.0e-6,
        xmax=0.032001,
        ymin=-1.0e-6,
        ymax=0.032001,
        zmin=-1.0e-6,
        zmax=0.032001,
    )

    remove_file(file_name)
    problem = UncollidedProblem(
        mesh=grid,
        num_groups=1,
        groupsets=[{"groups_from_to": [0, 0]}],
        xs_map=[{"block_ids": [0], "xs": xs}],
        point_sources=[point_source],
        near_source=[whole_domain],
        file_name=file_name,
        scattering_order=0,
    )

    scalar_flux = problem.GetScalarFluxFieldFunction()[0]
    minimum_scalar = volume_minimum(
        scalar_flux,
        whole_domain,
        FieldFunctionInterpolationVolume,
    )
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
    return point_errors, minimum_scalar


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
    point_errors = [result[0] for result in results]
    errors = [
        sum(error * error for error in mesh_errors) ** 0.5 / len(mesh_errors) ** 0.5
        for mesh_errors in point_errors
    ]
    max_errors = [max(mesh_errors) for mesh_errors in point_errors]
    minimum_scalar = min(result[1] for result in results)

    if rank == 0:
        print(f"Uncollided3DCoarseError={errors[0]:.8e}")
        print(f"Uncollided3DMediumError={errors[1]:.8e}")
        print(f"Uncollided3DFineError={errors[2]:.8e}")
        print(f"Uncollided3DFineMaxError={max_errors[2]:.8e}")
        print(f"Uncollided3DMinimumScalarFlux={minimum_scalar:.8e}")

    if not errors[1] < 0.8 * errors[0]:
        raise RuntimeError(f"3D coarse-to-medium convergence failed: {errors}")
    if not errors[2] < 0.8 * errors[1]:
        raise RuntimeError(f"3D medium-to-fine convergence failed: {errors}")
    if max_errors[2] > 0.05:
        raise RuntimeError(f"3D fine-mesh scalar-flux error is too large: {max_errors[2]}")
    if minimum_scalar < -1.0e-14:
        raise RuntimeError(f"Negative scalar flux remains: {minimum_scalar}")
