#!/usr/bin/env python3
"""Three-dimensional volumetric-source quadrature test."""

import math
import os
import sys

import numpy as np

if "opensn_console" not in globals():
    from mpi4py import MPI

    size = MPI.COMM_WORLD.size
    rank = MPI.COMM_WORLD.rank
    sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), "../../../../../")))
    from pyopensn.fieldfunc import FieldFunctionInterpolationPoint
    from pyopensn.logvol import RPPLogicalVolume
    from pyopensn.math import Vector3
    from pyopensn.mesh import OrthogonalMeshGenerator
    from pyopensn.solver import UncollidedProblem
    from pyopensn.source import VolumetricSource
    from pyopensn.xs import MultiGroupXS

from uncollided_unstructured_utils import point_value, relative_error, remove_file


if __name__ == "__main__":
    if size != 1:
        sys.exit(f"Expected one process, got {size}.")

    nodes = [-0.3, -0.1, 0.1, 0.3]
    grid = OrthogonalMeshGenerator(node_sets=[nodes, nodes, nodes]).Execute()
    grid.SetUniformBlockID(0)

    sigma_t = 0.7
    xs = MultiGroupXS()
    xs.CreateSimpleOneGroup(sigma_t=sigma_t, c=0.0)
    source_volume = RPPLogicalVolume(
        xmin=-0.1,
        xmax=0.1,
        ymin=-0.1,
        ymax=0.1,
        zmin=-0.1,
        zmax=0.1,
    )
    source = VolumetricSource(
        logical_volume=source_volume,
        group_strength=[1.0],
    )
    whole_domain = RPPLogicalVolume(
        xmin=-0.31,
        xmax=0.31,
        ymin=-0.31,
        ymax=0.31,
        zmin=-0.31,
        zmax=0.31,
    )

    file_name = "uncollided_volumetric_3d.h5"
    remove_file(file_name)
    problem = UncollidedProblem(
        mesh=grid,
        num_groups=1,
        groupsets=[{"groups_from_to": [0, 0]}],
        xs_map=[{"block_ids": [0], "xs": xs}],
        volumetric_sources=[source],
        volumetric_near_source=whole_domain,
        file_name=file_name,
        scattering_order=0,
    )

    sample = (0.25, 0.18, 0.22)
    value = point_value(
        problem.GetScalarFluxFieldFunction()[0],
        sample,
        FieldFunctionInterpolationPoint,
        Vector3,
    )
    reference = 0.0
    abscissae, weights = np.polynomial.legendre.leggauss(20)
    for xi, weight_x in zip(abscissae, weights):
        for eta, weight_y in zip(abscissae, weights):
            for zeta, weight_z in zip(abscissae, weights):
                point = (0.1 * xi, 0.1 * eta, 0.1 * zeta)
                distance = math.dist(point, sample)
                reference += (
                    0.1
                    * weight_x
                    * 0.1
                    * weight_y
                    * 0.1
                    * weight_z
                    * math.exp(-sigma_t * distance)
                    / (4.0 * math.pi * distance * distance)
                )
    error = relative_error(value, reference)

    if rank == 0:
        print(f"UncollidedVolumetric3DValue={value:.8e}")
        print(f"UncollidedVolumetric3DReference={reference:.8e}")
        print(f"UncollidedVolumetric3DQuadratureError={error:.8e}")

    remove_file(file_name)
    if error > 0.15:
        raise RuntimeError(f"3D volumetric-source quadrature error is too large: {error}")
