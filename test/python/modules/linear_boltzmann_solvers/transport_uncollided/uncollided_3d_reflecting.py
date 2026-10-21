#!/usr/bin/env python3
"""Analytic uncollided point-source test with three reflecting symmetry planes.

The neighbors of the source cell use exact_near_source_leakage.
"""

import importlib
import itertools
import math
import os
import sys

if "opensn_console" not in globals():
    from mpi4py import MPI

    size = MPI.COMM_WORLD.size
    rank = MPI.COMM_WORLD.rank
    sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), "../../../../../")))
    from pyopensn.aquad import GLCProductQuadrature3DXYZ
    from pyopensn.fieldfunc import FieldFunctionInterpolationPoint
    from pyopensn.logvol import RPPLogicalVolume
    from pyopensn.math import Vector3
    from pyopensn.mesh import OrthogonalMeshGenerator
    from pyopensn.solver import (
        DiscreteOrdinatesProblem,
        SteadyStateSourceSolver,
        UncollidedProblem,
        UncollidedSolver,
    )
    from pyopensn.source import PointSource
    from pyopensn.xs import MultiGroupXS

sys.path.append(os.path.dirname(__file__))
uncollided_utils = importlib.import_module("uncollided_unstructured_utils")
point_value = uncollided_utils.point_value
relative_error = uncollided_utils.relative_error
remove_file = uncollided_utils.remove_file


def analytic_flux(source, images, point, sigma_t):
    flux = 0.0
    for location in [source, *images]:
        distance = math.dist(point, location)
        flux += math.exp(-sigma_t * distance) / (4.0 * math.pi * distance**2)
    return flux


if __name__ == "__main__":
    if size != 1:
        sys.exit(f"Expected one process, got {size}.")

    nodes = [-1.0 + 0.05 * i for i in range(21)]
    grid = OrthogonalMeshGenerator(node_sets=[nodes, nodes, nodes]).Execute()
    grid.SetUniformBlockID(0)
    grid.SetOrthogonalBoundaries()

    sigma_t = 0.7
    xs = MultiGroupXS()
    xs.CreateSimpleOneGroup(sigma_t=sigma_t, c=0.0)

    # Reflections through the x = -1, y = -1, and z = -1 planes
    source = (-0.3589, -0.3865, -0.4421)
    images = [
        tuple(-2.0 - c if mirror else c for c, mirror in zip(source, mask))
        for mask in itertools.product((False, True), repeat=3)
        if any(mask)
    ]
    point_source = PointSource(location=list(source), strength=[1.0])
    # A small region around the point source
    near_source_region = RPPLogicalVolume(
        xmin=source[0] - 0.08,
        xmax=source[0] + 0.08,
        ymin=source[1] - 0.08,
        ymax=source[1] + 0.08,
        zmin=source[2] - 0.08,
        zmax=source[2] + 0.08,
    )
    boundary_conditions = [
        {"name": "xmin", "type": "reflecting"},
        {"name": "ymin", "type": "reflecting"},
        {"name": "zmin", "type": "reflecting"},
    ]

    file_name = "uncollided_3d_reflecting.h5"
    remove_file(file_name)
    uncollided = UncollidedProblem(
        mesh=grid,
        num_groups=1,
        groupsets=[{"groups_from_to": [0, 0]}],
        xs_map=[{"block_ids": [0], "xs": xs}],
        point_sources=[point_source],
        near_source=[near_source_region],
        exact_near_source_leakage=True,
        boundary_conditions=boundary_conditions,
        scattering_order=0,
    )
    uncollided_solver = UncollidedSolver(problem=uncollided, file_name=file_name)
    uncollided_solver.Initialize()
    uncollided_solver.Execute()

    sample_points = [
        (-0.74, -0.68, -0.83),
        (-0.21, -0.31, -0.12),
        (-0.62, -0.13, -0.71),
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

    quadrature = GLCProductQuadrature3DXYZ(
        n_polar=4,
        n_azimuthal=16,
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
        print(f"UncollidedReflecting3DMaxAnalyticError={max_analytic_error:.8e}")
        print(f"UncollidedReflecting3DRecombinationError={max_recombination_error:.8e}")
        print(f"UncollidedReflecting3DCombinedBalance={balance['balance']:.8e}")
        print(f"UncollidedReflecting3DMismatchRejected={int(mismatch_rejected)}")
        sys.stdout.flush()

    remove_file(file_name)
    if max_analytic_error > 0.02:
        raise RuntimeError(f"Reflecting analytic error is too large: {max_analytic_error}")
    if max_recombination_error > 1.0e-12:
        raise RuntimeError(
            f"Reflecting zero-scatter recombination changed the flux: {max_recombination_error}"
        )
    if not mismatch_rejected:
        raise RuntimeError("Mismatched reflecting boundaries were not rejected")
    # Directly projected images do not receive the paper's near-source correction.
    if abs(balance["balance"]) > 5.0e-3:
        raise RuntimeError(f"Reflecting combined balance failed: {balance['balance']}")
