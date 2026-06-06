#!/usr/bin/env python3
"""Cell-centroid volumetric-source quadrature and recombination test."""

import math
import os
import sys

import numpy as np

if "opensn_console" not in globals():
    from mpi4py import MPI

    size = MPI.COMM_WORLD.size
    rank = MPI.COMM_WORLD.rank
    sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), "../../../../../")))
    from pyopensn.aquad import GLCProductQuadrature2DXY
    from pyopensn.fieldfunc import (
        FieldFunctionInterpolationPoint,
        FieldFunctionInterpolationVolume,
    )
    from pyopensn.logvol import RPPLogicalVolume
    from pyopensn.math import Vector3, VectorSpatialFunction
    from pyopensn.mesh import OrthogonalMeshGenerator
    from pyopensn.solver import (
        DiscreteOrdinatesProblem,
        SteadyStateSourceSolver,
        UncollidedProblem,
    )
    from pyopensn.source import VolumetricSource
    from pyopensn.xs import MultiGroupXS

from uncollided_unstructured_utils import (
    point_value,
    relative_error,
    remove_file,
    volume_integral,
)


if __name__ == "__main__":
    if size != 1:
        sys.exit(f"Expected one process, got {size}.")

    spacing = 0.1
    nodes = [-1.0 + spacing * index for index in range(21)]
    grid = OrthogonalMeshGenerator(node_sets=[nodes, nodes]).Execute()
    grid.SetUniformBlockID(0)

    sigma_t = 0.7
    xs = MultiGroupXS()
    xs.CreateSimpleOneGroup(sigma_t=sigma_t, c=0.0)
    source_volume = RPPLogicalVolume(
        xmin=-0.2,
        xmax=0.2,
        ymin=-0.2,
        ymax=0.2,
        infz=True,
    )
    whole_domain = RPPLogicalVolume(
        xmin=-1.01,
        xmax=1.01,
        ymin=-1.01,
        ymax=1.01,
        infz=True,
    )

    def source_profile(point, num_groups):
        return [1.0 + point.x + 0.5 * point.y] * num_groups

    source = VolumetricSource(
        logical_volume=source_volume,
        func=VectorSpatialFunction(source_profile),
    )

    file_name = "uncollided_volumetric_2d.h5"
    remove_file(file_name)
    uncollided = UncollidedProblem(
        mesh=grid,
        num_groups=1,
        groupsets=[{"groups_from_to": [0, 0]}],
        xs_map=[{"block_ids": [0], "xs": xs}],
        volumetric_sources=[source],
        volumetric_near_source=whole_domain,
        file_name=file_name,
        scattering_order=0,
    )

    sample = (0.72, 0.13, 0.0)
    value = point_value(
        uncollided.GetScalarFluxFieldFunction()[0],
        sample,
        FieldFunctionInterpolationPoint,
        Vector3,
    )
    reference = 0.0
    abscissae, weights = np.polynomial.legendre.leggauss(32)
    for xi, weight_x in zip(abscissae, weights):
        for eta, weight_y in zip(abscissae, weights):
            x = 0.2 * xi
            y = 0.2 * eta
            distance = math.dist((x, y), sample[:2])
            source_value = 1.0 + x + 0.5 * y
            reference += (
                0.2
                * weight_x
                * 0.2
                * weight_y
                * source_value
                * math.exp(-sigma_t * distance)
                / (2.0 * math.pi * distance)
            )
    quadrature_error = relative_error(value, reference)

    quadrature = GLCProductQuadrature2DXY(
        n_polar=2,
        n_azimuthal=12,
        scattering_order=0,
    )
    problem = DiscreteOrdinatesProblem(
        mesh=grid,
        num_groups=1,
        groupsets=[
            {
                "groups_from_to": [0, 0],
                "angular_quadrature": quadrature,
            }
        ],
        xs_map=[{"block_ids": [0], "xs": xs}],
        uncollided_flux=file_name,
    )
    solver = SteadyStateSourceSolver(problem=problem)
    solver.Initialize()
    solver.Execute()
    recombined = point_value(
        problem.GetScalarFluxFieldFunction()[0],
        sample,
        FieldFunctionInterpolationPoint,
        Vector3,
    )
    recombination_error = relative_error(recombined, value)

    scattering_xs = MultiGroupXS()
    scattering_xs.CreateSimpleOneGroup(sigma_t=sigma_t, c=0.55)
    collided_problem = DiscreteOrdinatesProblem(
        mesh=grid,
        num_groups=1,
        groupsets=[
            {
                "groups_from_to": [0, 0],
                "angular_quadrature": quadrature,
                "inner_linear_method": "petsc_gmres",
                "l_abs_tol": 1.0e-10,
                "l_max_its": 200,
            }
        ],
        xs_map=[{"block_ids": [0], "xs": scattering_xs}],
        uncollided_flux=file_name,
    )
    collided_solver = SteadyStateSourceSolver(problem=collided_problem, compute_balance=True)
    collided_solver.Initialize()
    collided_solver.Execute()
    uncollided_integral = volume_integral(
        uncollided.GetScalarFluxFieldFunction()[0],
        whole_domain,
        FieldFunctionInterpolationVolume,
    )
    total_integral = volume_integral(
        collided_problem.GetScalarFluxFieldFunction()[0],
        whole_domain,
        FieldFunctionInterpolationVolume,
    )
    collided_fraction = (total_integral - uncollided_integral) / total_integral
    combined_balance = collided_solver.ComputeBalanceTable()["balance"]

    if rank == 0:
        print(f"UncollidedVolumetric2DValue={value:.8e}")
        print(f"UncollidedVolumetric2DReference={reference:.8e}")
        print(f"UncollidedVolumetric2DQuadratureError={quadrature_error:.8e}")
        print(f"UncollidedVolumetric2DRecombinationError={recombination_error:.8e}")
        print(f"UncollidedVolumetric2DCombinedBalance={combined_balance:.8e}")
        print(f"UncollidedVolumetric2DCollidedFraction={collided_fraction:.8e}")

    remove_file(file_name)
    if quadrature_error > 0.08:
        raise RuntimeError(f"Volumetric-source quadrature error is too large: {quadrature_error}")
    if recombination_error > 1.0e-12:
        raise RuntimeError(f"Volumetric-source recombination failed: {recombination_error}")
    if abs(combined_balance) > 2.0e-5:
        raise RuntimeError(f"Volumetric-source combined balance failed: {combined_balance}")
    if collided_fraction < 0.1:
        raise RuntimeError(
            f"Volumetric first-collision source is too small: {collided_fraction}"
        )
