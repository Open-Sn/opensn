#!/usr/bin/env python3
"""Exact source-cell average of the uncollided flux on an orthogonal 3D mesh.

The flux in the cell containing a point source is represented by its exact cell
average. In a homogeneous cell, integrating the uncollided flux along rays from
the source gives

    int_V phi dV = S / (4 pi) * int (1 - exp(-sigma_t R)) / sigma_t dOmega,

where R is the distance from the source to the cell boundary in direction
Omega. This test evaluates the reference independently, as Cartesian integrals
over each face using dOmega = h / r^3 dA, where h is the distance from the
source to the face plane. The composite Gauss-Legendre rules converge to
round-off. The source sits close to a face, where the earlier quadrature
projection was least accurate.

A source on a vacuum boundary face emits only the half of its strength that
enters the domain into the cell, so its cell average is half that of the
corresponding interior source cut at the face. The face containing the source
subtends no angle and drops out of the reference.
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
relative_error = uncollided_utils.relative_error
remove_file = uncollided_utils.remove_file
box_angular_integral = uncollided_utils.box_angular_integral

SIGMA_T = 2.0
STRENGTH = 1.0


def reference_average(bounds, source):
    """Exact cell-average flux of a point source in an axis-aligned box."""
    volume_angle_integral = box_angular_integral(
        lambda r: -math.expm1(-SIGMA_T * r) / SIGMA_T, source, bounds
    )
    cell_volume = math.prod(b[1] - b[0] for b in bounds)
    return STRENGTH / (4.0 * math.pi) * volume_angle_integral / cell_volume


def computed_average(source, cell_bounds, file_name):
    nodes = [0.0, 1.0, 2.0, 3.0]
    grid = OrthogonalMeshGenerator(node_sets=[nodes] * 3).Execute()
    grid.SetUniformBlockID(0)

    xs = MultiGroupXS()
    xs.CreateSimpleOneGroup(sigma_t=SIGMA_T, c=0.0)

    # Ray-trace the whole mesh so the result does not depend on the bulk sweep, which is
    # inaccurate on these optically thick (sigma_t h = 2) cells.
    region = RPPLogicalVolume(xmin=-1.0, xmax=4.0, ymin=-1.0, ymax=4.0, zmin=-1.0, zmax=4.0)

    remove_file(file_name)
    problem = UncollidedProblem(
        mesh=grid,
        num_groups=1,
        groupsets=[{"groups_from_to": [0, 0]}],
        xs_map=[{"block_ids": [0], "xs": xs}],
        point_sources=[PointSource(location=list(source), strength=[STRENGTH])],
        near_source=[region],
        scattering_order=0,
    )
    solver = UncollidedSolver(problem=problem, file_name=file_name)
    solver.Initialize()
    solver.Execute()

    # The source-cell flux is constant, so any interior point gives its average.
    centroid = [0.5 * (b[0] + b[1]) for b in cell_bounds]
    scalar_flux = problem.GetScalarFluxFieldFunction()[0]
    value = point_value(scalar_flux, centroid, FieldFunctionInterpolationPoint, Vector3)
    remove_file(file_name)
    return value


if __name__ == "__main__":
    if size != 1:
        sys.exit(f"Expected one process, got {size}.")

    cases = {
        "Interior": ((1.37, 1.62, 1.06), [(1.0, 2.0)] * 3),
        "Boundary": ((1.37, 1.62, 0.0), [(1.0, 2.0), (1.0, 2.0), (0.0, 1.0)]),
    }
    errors = {}
    for name, (source, cell_bounds) in cases.items():
        value = computed_average(source, cell_bounds, f"uncollided_source_cell_{name}.h5")
        reference = reference_average(cell_bounds, source)
        errors[name] = relative_error(value, reference)

    if rank == 0:
        for name, error in errors.items():
            print(f"UncollidedSourceCell{name}Error={error:.8e}")
