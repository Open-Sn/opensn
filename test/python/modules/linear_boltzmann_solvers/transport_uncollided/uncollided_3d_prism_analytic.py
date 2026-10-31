#!/usr/bin/env python3
"""Homogeneous 3D prism-mesh uncollided flux and P1 moment test.

The mesh extrudes an unstructured triangle mesh, so cells have triangular and
quadrilateral faces. Scalar flux and P1 moments are compared with the exact
point-source solution, and the global outflow with the exact uncollided escape
fraction of the box. Near the singularity individual PWLD coefficients of the
ray-traced cells may be negative (see the theory manual), so positivity is
checked in the swept bulk region outside the near-source box.
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
    from pyopensn.aquad import GLCProductQuadrature3DXYZ
    from pyopensn.logvol import RPPLogicalVolume
    from pyopensn.fieldfunc import (
        FieldFunctionInterpolationPoint,
        FieldFunctionInterpolationVolume,
    )
    from pyopensn.math import Vector3
    from pyopensn.mesh import ExtruderMeshGenerator, FromFileMeshGenerator
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
extruded_mesh = uncollided_utils.extruded_mesh
point_value = uncollided_utils.point_value
relative_error = uncollided_utils.relative_error
remove_file = uncollided_utils.remove_file
uncollided_3d = uncollided_utils.uncollided_3d
uncollided_escape_3d_box = uncollided_utils.uncollided_escape_3d_box
volume_integral = uncollided_utils.volume_integral
volume_minimum = uncollided_utils.volume_minimum


if __name__ == "__main__":
    if size != 1:
        sys.exit(f"Expected one process, got {size}.")

    grid = extruded_mesh(
        ExtruderMeshGenerator, FromFileMeshGenerator, "triangle_mesh2x2.obj", num_layers=16
    )
    grid.SetUniformBlockID(0)

    sigma_t = 0.7
    xs = MultiGroupXS()
    xs.CreateSimpleOneGroup(sigma_t=sigma_t, c=0.0)

    source = (0.0441723, -0.0330477, 1.11)
    strength = 1.0
    point_source = PointSource(location=list(source), strength=[strength])
    whole_domain = RPPLogicalVolume(
        xmin=-1.01,
        xmax=1.01,
        ymin=-1.01,
        ymax=1.01,
        zmin=-0.01,
        zmax=2.01,
    )
    # Small region around the source
    half_width = 0.15
    near_source_region = RPPLogicalVolume(
        xmin=source[0] - half_width,
        xmax=source[0] + half_width,
        ymin=source[1] - half_width,
        ymax=source[1] + half_width,
        zmin=source[2] - half_width,
        zmax=source[2] + half_width,
    )
    # Six slabs that together cover the domain outside the near-source box
    domain = [(-1.01, 1.01), (-1.01, 1.01), (-0.01, 2.01)]
    bulk_slabs = []
    for axis in range(3):
        for side in range(2):
            bounds = list(domain)
            if side == 0:
                bounds[axis] = (domain[axis][0], source[axis] - half_width)
            else:
                bounds[axis] = (source[axis] + half_width, domain[axis][1])
            bulk_slabs.append(
                RPPLogicalVolume(
                    xmin=bounds[0][0],
                    xmax=bounds[0][1],
                    ymin=bounds[1][0],
                    ymax=bounds[1][1],
                    zmin=bounds[2][0],
                    zmax=bounds[2][1],
                )
            )

    file_name = "uncollided_3d_prism_analytic.h5"
    remove_file(file_name)
    uncollided = UncollidedProblem(
        mesh=grid,
        num_groups=1,
        groupsets=[{"groups_from_to": [0, 0]}],
        xs_map=[{"block_ids": [0], "xs": xs}],
        point_sources=[point_source],
        near_source=[near_source_region],
        scattering_order=1,
    )
    uncollided_solver = UncollidedSolver(problem=uncollided, file_name=file_name)
    uncollided_solver.Initialize()
    uncollided_solver.Execute()

    uncollided_scalar = uncollided.GetScalarFluxFieldFunction()[0]
    minimum_scalar = min(
        volume_minimum(uncollided_scalar, slab, FieldFunctionInterpolationVolume)
        for slab in bulk_slabs
    )
    scalar_integral = volume_integral(
        uncollided_scalar,
        whole_domain,
        FieldFunctionInterpolationVolume,
    )
    computed_outflow = strength - sigma_t * scalar_integral
    reference_outflow = uncollided_escape_3d_box(
        sigma_t, source, [(-1.0, 1.0), (-1.0, 1.0), (0.0, 2.0)]
    )
    outflow_error = relative_error(computed_outflow, reference_outflow)
    quadrature = GLCProductQuadrature3DXYZ(
        n_polar=4,
        n_azimuthal=16,
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
        (0.46, 0.11, 1.37),
        (-0.39, 0.43, 0.82),
        (0.22, -0.58, 1.21),
    ]
    max_scalar_error = 0.0
    max_p1_error = 0.0
    max_recombination_error = 0.0
    for point in sample_points:
        scalar = point_value(moments[0], point, FieldFunctionInterpolationPoint, Vector3)
        original_scalar = point_value(
            uncollided_scalar, point, FieldFunctionInterpolationPoint, Vector3
        )
        exact_scalar = uncollided_3d(strength, sigma_t, source, point)
        max_scalar_error = max(max_scalar_error, relative_error(scalar, exact_scalar))
        max_recombination_error = max(
            max_recombination_error,
            abs(scalar - original_scalar) / max(abs(original_scalar), 1.0e-14),
        )

        dx, dy, dz = (p - q for p, q in zip(point, source))
        radius = math.sqrt(dx * dx + dy * dy + dz * dz)
        exact_p1 = (
            exact_scalar * dy / radius,
            exact_scalar * dz / radius,
            exact_scalar * dx / radius,
        )
        for moment, reference in zip(moments[1:], exact_p1):
            value = point_value(moment, point, FieldFunctionInterpolationPoint, Vector3)
            scale = max(abs(exact_scalar), 1.0e-14)
            max_p1_error = max(max_p1_error, abs(value - reference) / scale)

    if rank == 0:
        print(f"Uncollided3DPrismMaxScalarRelativeError={max_scalar_error:.8e}")
        print(f"Uncollided3DPrismMaxP1ScaledError={max_p1_error:.8e}")
        print(f"Uncollided3DPrismRecombinationError={max_recombination_error:.8e}")
        print(f"Uncollided3DPrismBulkMinimumScalarFlux={minimum_scalar:.8e}")
        print(f"Uncollided3DPrismOutflowRelativeError={outflow_error:.8e}")
        sys.stdout.flush()

    remove_file(file_name)
    if max_scalar_error > 0.05:
        raise RuntimeError(f"3D prism scalar-flux error is too large: {max_scalar_error}")
    if max_p1_error > 0.05:
        raise RuntimeError(f"3D prism P1-moment error is too large: {max_p1_error}")
    if max_recombination_error > 1.0e-12:
        raise RuntimeError(
            f"Zero-scatter recombination changed scalar flux: {max_recombination_error}"
        )
    if minimum_scalar < -1.0e-14:
        raise RuntimeError(f"Negative bulk-region scalar flux remains: {minimum_scalar}")
    if outflow_error > 2.0e-3:
        raise RuntimeError(f"3D prism global outflow error is too large: {outflow_error}")
