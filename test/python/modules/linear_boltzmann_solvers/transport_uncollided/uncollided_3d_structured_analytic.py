#!/usr/bin/env python3
"""Verify 3D structured-mesh cell-average positivity and global outflow.

The outflow is compared with the exact uncollided escape fraction of the box. Near the
point-source singularity individual PWLD coefficients may be negative (see the theory
manual), so positivity is checked on cell averages. On these orthogonal hexahedra the
trilinear value at a cell centroid equals the cell average.
"""

import importlib
import os
import sys

if "opensn_console" not in globals():
    from mpi4py import MPI

    rank = MPI.COMM_WORLD.rank
    size = MPI.COMM_WORLD.size
    sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), "../../../../../")))
    from pyopensn.fieldfunc import (
        FieldFunctionInterpolationPoint,
        FieldFunctionInterpolationVolume,
    )
    from pyopensn.logvol import RPPLogicalVolume
    from pyopensn.math import Vector3
    from pyopensn.mesh import OrthogonalMeshGenerator
    from pyopensn.solver import UncollidedProblem, UncollidedSolver
    from pyopensn.source import PointSource
    from pyopensn.xs import MultiGroupXS

sys.path.append(os.path.dirname(__file__))
uncollided_utils = importlib.import_module("uncollided_unstructured_utils")
remove_file = uncollided_utils.remove_file
uncollided_escape_3d_box = uncollided_utils.uncollided_escape_3d_box
point_value = uncollided_utils.point_value
volume_integral = uncollided_utils.volume_integral


if __name__ == "__main__":
    if size != 1:
        sys.exit(f"Expected one process, got {size}.")

    nodes = [-1.0, -0.5, 0.0, 0.5, 1.0]
    grid = OrthogonalMeshGenerator(node_sets=[nodes, nodes, nodes]).Execute()
    grid.SetUniformBlockID(0)

    xs = MultiGroupXS()
    xs.CreateSimpleOneGroup(sigma_t=0.8, c=0.0)

    whole_domain = RPPLogicalVolume(
        xmin=-1.01,
        xmax=1.01,
        ymin=-1.01,
        ymax=1.01,
        zmin=-1.01,
        zmax=1.01,
    )

    source = (0.14, -0.19, 0.23)
    file_name = "uncollided_3d_structured_analytic.h5"
    remove_file(file_name)
    problem = UncollidedProblem(
        mesh=grid,
        num_groups=1,
        groupsets=[{"groups_from_to": [0, 0]}],
        xs_map=[{"block_ids": [0], "xs": xs}],
        point_sources=[PointSource(location=list(source), strength=[1.0])],
        near_source=[whole_domain],
        scattering_order=0,
    )
    solver = UncollidedSolver(problem=problem, file_name=file_name)
    solver.Initialize()
    solver.Execute()

    scalar_flux = problem.GetScalarFluxFieldFunction()[0]
    centroids = [0.5 * (a + b) for a, b in zip(nodes[:-1], nodes[1:])]
    minimum_scalar = min(
        point_value(scalar_flux, (x, y, z), FieldFunctionInterpolationPoint, Vector3)
        for x in centroids
        for y in centroids
        for z in centroids
    )
    negative_scalar_violation = max(0.0, -minimum_scalar)
    scalar_integral = volume_integral(scalar_flux, whole_domain, FieldFunctionInterpolationVolume)
    computed_outflow = 1.0 - 0.8 * scalar_integral
    reference_outflow = uncollided_escape_3d_box(0.8, source, [(-1.0, 1.0)] * 3)
    outflow_error = abs(computed_outflow - reference_outflow) / abs(reference_outflow)
    remove_file(file_name)

    if rank == 0:
        print(f"Uncollided3DStructuredNegativeCellAverageViolation={negative_scalar_violation:.8e}")
        print(f"Uncollided3DStructuredOutflowRelativeError={outflow_error:.8e}")
        sys.stdout.flush()

    if minimum_scalar < -1.0e-14:
        raise RuntimeError(f"Negative cell-average scalar flux remains: {minimum_scalar}")
    if outflow_error > 2.0e-2:
        raise RuntimeError(f"3D structured outflow error is too large: {outflow_error}")
