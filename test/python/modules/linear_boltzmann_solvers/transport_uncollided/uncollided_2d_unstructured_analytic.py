#!/usr/bin/env python3
"""Homogeneous 2D unstructured-mesh uncollided flux and P1 moment test."""

import importlib
import math
import os
import sys

if "opensn_console" not in globals():
    from mpi4py import MPI

    size = MPI.COMM_WORLD.size
    rank = MPI.COMM_WORLD.rank
    sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), "../../../../../")))
    from pyopensn.aquad import GLCProductQuadrature2DXY
    from pyopensn.logvol import RPPLogicalVolume
    from pyopensn.fieldfunc import (
        FieldFunctionInterpolationPoint,
        FieldFunctionInterpolationVolume,
    )
    from pyopensn.math import Vector3
    from pyopensn.mesh import FromFileMeshGenerator
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
mesh_path = uncollided_utils.mesh_path
point_value = uncollided_utils.point_value
relative_error = uncollided_utils.relative_error
remove_file = uncollided_utils.remove_file
uncollided_2d = uncollided_utils.uncollided_2d
uncollided_escape_2d_rectangle = uncollided_utils.uncollided_escape_2d_rectangle
volume_integral = uncollided_utils.volume_integral
volume_minimum = uncollided_utils.volume_minimum


if __name__ == "__main__":
    if size != 1:
        sys.exit(f"Expected one process, got {size}.")

    grid = FromFileMeshGenerator(filename=mesh_path("triangle_mesh2x2_fine.obj")).Execute()
    grid.SetUniformBlockID(0)

    sigma_t = 0.7
    xs = MultiGroupXS()
    xs.CreateSimpleOneGroup(sigma_t=sigma_t, c=0.0)

    source = (0.037, -0.041, 0.0)
    strength = 1.0
    point_source = PointSource(location=list(source), strength=[strength])
    whole_domain = RPPLogicalVolume(
        xmin=-1.01,
        xmax=1.01,
        ymin=-1.01,
        ymax=1.01,
        infz=True,
    )

    file_name = "uncollided_2d_unstructured_analytic.h5"
    remove_file(file_name)
    uncollided = UncollidedProblem(
        mesh=grid,
        num_groups=1,
        groupsets=[{"groups_from_to": [0, 0]}],
        xs_map=[{"block_ids": [0], "xs": xs}],
        point_sources=[point_source],
        near_source=[whole_domain],
        scattering_order=1,
    )
    uncollided_solver = UncollidedSolver(problem=uncollided, file_name=file_name)
    uncollided_solver.Initialize()
    uncollided_solver.Execute()

    uncollided_scalar = uncollided.GetScalarFluxFieldFunction()[0]
    minimum_scalar = volume_minimum(
        uncollided_scalar,
        whole_domain,
        FieldFunctionInterpolationVolume,
    )
    scalar_integral = volume_integral(
        uncollided_scalar,
        whole_domain,
        FieldFunctionInterpolationVolume,
    )
    computed_outflow = strength - sigma_t * scalar_integral
    reference_outflow = uncollided_escape_2d_rectangle(
        sigma_t,
        source,
        (-1.0, 1.0, -1.0, 1.0),
    )
    outflow_error = relative_error(computed_outflow, reference_outflow)
    quadrature = GLCProductQuadrature2DXY(
        n_polar=2,
        n_azimuthal=12,
        scattering_order=1,
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
    moments = problem.GetScalarFluxFieldFunction(only_scalar_flux=False)[0]

    sample_points = [
        (0.46, 0.11, 0.0),
        (-0.39, 0.43, 0.0),
        (0.22, -0.58, 0.0),
    ]
    max_scalar_error = 0.0
    max_p1_error = 0.0
    max_recombination_error = 0.0
    for point in sample_points:
        scalar = point_value(moments[0], point, FieldFunctionInterpolationPoint, Vector3)
        original_scalar = point_value(
            uncollided_scalar, point, FieldFunctionInterpolationPoint, Vector3
        )
        exact_scalar = uncollided_2d(strength, sigma_t, source, point)
        max_scalar_error = max(max_scalar_error, relative_error(scalar, exact_scalar))
        max_recombination_error = max(
            max_recombination_error,
            abs(scalar - original_scalar) / max(abs(original_scalar), 1.0e-14),
        )

        dx = point[0] - source[0]
        dy = point[1] - source[1]
        radius = math.hypot(dx, dy)
        exact_p1 = (exact_scalar * dy / radius, exact_scalar * dx / radius)
        for moment, reference in zip(moments[1:], exact_p1):
            value = point_value(moment, point, FieldFunctionInterpolationPoint, Vector3)
            scale = max(abs(exact_scalar), 1.0e-14)
            max_p1_error = max(max_p1_error, abs(value - reference) / scale)

    if rank == 0:
        print(f"Uncollided2DMaxScalarRelativeError={max_scalar_error:.8e}")
        print(f"Uncollided2DMaxP1ScaledError={max_p1_error:.8e}")
        print(f"Uncollided2DRecombinationError={max_recombination_error:.8e}")
        print(f"Uncollided2DMinimumScalarFlux={minimum_scalar:.8e}")
        print(f"Uncollided2DOutflowRelativeError={outflow_error:.8e}")
        sys.stdout.flush()

    remove_file(file_name)
    if max_scalar_error > 5.0e-4:
        raise RuntimeError(f"2D scalar-flux error is too large: {max_scalar_error}")
    if max_p1_error > 1.0e-3:
        raise RuntimeError(f"2D P1-moment error is too large: {max_p1_error}")
    if max_recombination_error > 1.0e-12:
        raise RuntimeError(
            f"Zero-scatter recombination changed scalar flux: {max_recombination_error}"
        )
    if minimum_scalar < -1.0e-14:
        raise RuntimeError(f"Negative scalar flux remains: {minimum_scalar}")
    if outflow_error > 1.0e-3:
        raise RuntimeError(f"2D global outflow error is too large: {outflow_error}")
