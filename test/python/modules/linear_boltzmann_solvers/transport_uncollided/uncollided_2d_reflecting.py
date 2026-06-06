#!/usr/bin/env python3
"""Analytic uncollided point-source test with two reflecting symmetry planes."""

import math
import os
import sys

if "opensn_console" not in globals():
    from mpi4py import MPI

    size = MPI.COMM_WORLD.size
    rank = MPI.COMM_WORLD.rank
    sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), "../../../../../")))
    from pyopensn.aquad import GLCProductQuadrature2DXY
    from pyopensn.fieldfunc import FieldFunctionInterpolationPoint
    from pyopensn.logvol import RPPLogicalVolume
    from pyopensn.math import Vector3
    from pyopensn.mesh import FromFileMeshGenerator
    from pyopensn.solver import (
        DiscreteOrdinatesProblem,
        SteadyStateSourceSolver,
        UncollidedProblem,
    )
    from pyopensn.source import PointSource
    from pyopensn.xs import MultiGroupXS

from uncollided_unstructured_utils import mesh_path, point_value, relative_error, remove_file


def analytic_flux(source, images, point, sigma_t):
    flux = 0.0
    for location in [source, *images]:
        distance = math.hypot(point[0] - location[0], point[1] - location[1])
        flux += math.exp(-sigma_t * distance) / (2.0 * math.pi * distance)
    return flux


if __name__ == "__main__":
    if size != 1:
        sys.exit(f"Expected one process, got {size}.")

    grid = FromFileMeshGenerator(filename=mesh_path("triangle_mesh2x2_fine.obj")).Execute()
    grid.SetUniformBlockID(0)
    grid.SetOrthogonalBoundaries()

    sigma_t = 0.7
    xs = MultiGroupXS()
    xs.CreateSimpleOneGroup(sigma_t=sigma_t, c=0.0)

    source = (-0.37, -0.41, 0.0)
    images = [
        (-1.63, -0.41, 0.0),
        (-0.37, -1.59, 0.0),
        (-1.63, -1.59, 0.0),
    ]
    point_source = PointSource(location=list(source), strength=[1.0])
    whole_domain = RPPLogicalVolume(
        xmin=-1.01,
        xmax=1.01,
        ymin=-1.01,
        ymax=1.01,
        infz=True,
    )
    boundary_conditions = [
        {"name": "xmin", "type": "reflecting"},
        {"name": "ymin", "type": "reflecting"},
    ]

    file_name = "uncollided_2d_reflecting.h5"
    remove_file(file_name)
    uncollided = UncollidedProblem(
        mesh=grid,
        num_groups=1,
        groupsets=[{"groups_from_to": [0, 0]}],
        xs_map=[{"block_ids": [0], "xs": xs}],
        point_sources=[point_source],
        near_source=[whole_domain],
        boundary_conditions=boundary_conditions,
        file_name=file_name,
        scattering_order=0,
    )

    sample_points = [
        (-0.74, -0.68, 0.0),
        (0.21, -0.31, 0.0),
        (-0.12, 0.46, 0.0),
    ]
    scalar_flux = uncollided.GetScalarFluxFieldFunction()[0]
    max_analytic_error = 0.0
    for point in sample_points:
        value = point_value(
            scalar_flux,
            point,
            FieldFunctionInterpolationPoint,
            Vector3,
        )
        reference = analytic_flux(source, images, point, sigma_t)
        max_analytic_error = max(max_analytic_error, relative_error(value, reference))

    projected_file_name = "uncollided_2d_reflecting_projected.h5"
    remove_file(projected_file_name)
    projected_uncollided = UncollidedProblem(
        mesh=grid,
        num_groups=1,
        groupsets=[{"groups_from_to": [0, 0]}],
        xs_map=[{"block_ids": [0], "xs": xs}],
        point_sources=[point_source],
        near_source=[whole_domain],
        boundary_conditions=boundary_conditions,
        file_name=projected_file_name,
        scattering_order=0,
        project_reflected_image_sources=True,
        reflected_image_projection_threads=2,
    )
    projected_scalar_flux = projected_uncollided.GetScalarFluxFieldFunction()[0]
    max_projected_analytic_error = 0.0
    for point in sample_points:
        value = point_value(
            projected_scalar_flux,
            point,
            FieldFunctionInterpolationPoint,
            Vector3,
        )
        reference = analytic_flux(source, images, point, sigma_t)
        max_projected_analytic_error = max(
            max_projected_analytic_error, relative_error(value, reference)
        )

    quadrature = GLCProductQuadrature2DXY(
        n_polar=2,
        n_azimuthal=24,
        scattering_order=0,
    )
    mismatch_rejected = False
    try:
        DiscreteOrdinatesProblem(
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
    except ValueError:
        mismatch_rejected = True

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
        boundary_conditions=boundary_conditions,
        uncollided_flux=file_name,
    )
    solver = SteadyStateSourceSolver(problem=problem, compute_balance=True)
    solver.Initialize()
    solver.Execute()

    recombined_flux = problem.GetScalarFluxFieldFunction()[0]
    max_recombination_error = 0.0
    for point in sample_points:
        original = point_value(
            scalar_flux,
            point,
            FieldFunctionInterpolationPoint,
            Vector3,
        )
        recombined = point_value(
            recombined_flux,
            point,
            FieldFunctionInterpolationPoint,
            Vector3,
        )
        max_recombination_error = max(
            max_recombination_error,
            abs(recombined - original) / max(abs(original), 1.0e-14),
        )

    balance = solver.ComputeBalanceTable()
    if rank == 0:
        print(f"UncollidedReflecting2DMaxAnalyticError={max_analytic_error:.8e}")
        print(
            "UncollidedReflecting2DProjectedMaxAnalyticError="
            f"{max_projected_analytic_error:.8e}"
        )
        print(f"UncollidedReflecting2DRecombinationError={max_recombination_error:.8e}")
        print(f"UncollidedReflecting2DCombinedBalance={balance['balance']:.8e}")
        print(f"UncollidedReflecting2DMismatchRejected={int(mismatch_rejected)}")
        sys.stdout.flush()

    remove_file(file_name)
    remove_file(projected_file_name)
    if max_analytic_error > 2.0e-3:
        raise RuntimeError(f"Reflecting analytic error is too large: {max_analytic_error}")
    if max_projected_analytic_error > 2.0e-3:
        raise RuntimeError(
            "Projected reflecting analytic error is too large: "
            f"{max_projected_analytic_error}"
        )
    if max_recombination_error > 1.0e-12:
        raise RuntimeError(
            f"Reflecting zero-scatter recombination changed the flux: {max_recombination_error}"
        )
    if not mismatch_rejected:
        raise RuntimeError("Mismatched reflecting boundaries were not rejected")
    if abs(balance["balance"]) > 2.0e-5:
        raise RuntimeError(f"Reflecting combined balance failed: {balance['balance']}")
