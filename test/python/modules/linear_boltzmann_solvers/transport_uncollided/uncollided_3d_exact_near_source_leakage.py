#!/usr/bin/env python3
"""Exact near-source leakage for a point source close to a vertex of its cell.

When a point source lies close to an edge or vertex of its cell, the neighboring
cells are close to the singular flux too, and their spatial-quadrature leakages
carry a large error into every downstream cell. With exact_near_source_leakage
those cells are integrated exactly over direction. The scalar flux at 40 points
away from the source is compared with the exact point-source solution in a
homogeneous medium. With the option the error matches that of a source at the
cell center on this mesh (about 2.5% RMS), which is the bulk-sweep
discretization error. Without it the error is about 23% RMS; that value is
printed for comparison but not checked.
"""

import importlib
import math
import os
import random
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
remove_file = uncollided_utils.remove_file
uncollided_3d = uncollided_utils.uncollided_3d

SIGMA_T = 0.25
H = 0.05
LOWER = (-0.7, 0.0, -0.4)
UPPER = (0.7, 0.3, 0.2)
# 0.03 cell widths from the vertex (-0.6, 0.1, 0.0) of its cell
SOURCE = (-0.5985, 0.1015, 0.0015)


def sample_points():
    generator = random.Random(1)
    points = []
    while len(points) < 40:
        point = tuple(generator.uniform(LOWER[i] + 0.01, UPPER[i] - 0.01) for i in range(3))
        if math.dist(point, SOURCE) > 0.3:
            points.append(point)
    return points


def relative_errors(grid, exact_near_source_leakage, points, file_name):
    xs = MultiGroupXS()
    xs.CreateSimpleOneGroup(sigma_t=SIGMA_T, c=0.0)
    near_source_region = RPPLogicalVolume(
        xmin=SOURCE[0] - 0.08,
        xmax=SOURCE[0] + 0.08,
        ymin=SOURCE[1] - 0.08,
        ymax=SOURCE[1] + 0.08,
        zmin=SOURCE[2] - 0.08,
        zmax=SOURCE[2] + 0.08,
    )

    remove_file(file_name)
    problem = UncollidedProblem(
        mesh=grid,
        num_groups=1,
        groupsets=[{"groups_from_to": [0, 0]}],
        xs_map=[{"block_ids": [0], "xs": xs}],
        point_sources=[PointSource(location=list(SOURCE), strength=[1.0])],
        near_source=[near_source_region],
        scattering_order=0,
        exact_near_source_leakage=exact_near_source_leakage,
    )
    solver = UncollidedSolver(problem=problem, file_name=file_name)
    solver.Initialize()
    solver.Execute()

    scalar_flux = problem.GetScalarFluxFieldFunction()[0]
    errors = []
    for point in points:
        value = point_value(scalar_flux, point, FieldFunctionInterpolationPoint, Vector3)
        reference = uncollided_3d(1.0, SIGMA_T, SOURCE, point)
        errors.append((value - reference) / reference)
    remove_file(file_name)
    return errors


def rms(values):
    return math.sqrt(sum(value * value for value in values) / len(values))


if __name__ == "__main__":
    if size != 1:
        sys.exit(f"Expected one process, got {size}.")

    node_sets = [
        [LOWER[i] + H * k for k in range(round((UPPER[i] - LOWER[i]) / H) + 1)] for i in range(3)
    ]
    grid = OrthogonalMeshGenerator(node_sets=node_sets).Execute()
    grid.SetUniformBlockID(0)

    points = sample_points()
    exact = relative_errors(grid, True, points, "uncollided_exact_leakage_on.h5")
    quadrature = relative_errors(grid, False, points, "uncollided_exact_leakage_off.h5")

    exact_rms = rms(exact)
    exact_max = max(abs(error) for error in exact)
    if rank == 0:
        print(f"UncollidedExactLeakageRMSError={exact_rms:.8e}")
        print(f"UncollidedExactLeakageMaxError={exact_max:.8e}")
        print(f"UncollidedQuadratureLeakageRMSError={rms(quadrature):.8e}")
        sys.stdout.flush()

    if exact_rms > 0.03:
        raise RuntimeError(f"Exact near-source leakage RMS error is too large: {exact_rms}")
    if exact_max > 0.06:
        raise RuntimeError(f"Exact near-source leakage max error is too large: {exact_max}")
