#!/usr/bin/env python3
"""Uncollided ray tracing on an orthogonal mesh, where rays pass exactly through edges and
vertices, across two materials.

The box [0, 2000] x [0, 2000] x [-100, 2000] (in cm, so that the coordinates are large compared
with the thinnest cells) has 100 x 100 columns of hexahedra. Below z = 0 it is a dense ground
layer (sigma_t = 0.15 / cm) on cells that grade from 40 cm thick to 1 cm thick at its surface;
above it is a thin medium (sigma_t = 1e-4 / cm). The point source at (1030, 1010, 100) lies
inside a cell. On this regular grid many rays from quadrature points to the source pass exactly
through cell edges and vertices. The earlier ray tracer then nudged rays, lost them, and
stopped with "could not locate the cell ahead of a ray". The whole box is ray traced, with the
exact near-source leakage.

The outflow, one minus the absorption in the box, is compared with the exact escape fraction:
the integral over the box boundary of exp(-tau) h / (4 pi r^3), where tau is the optical depth
of the straight line from the source to the boundary point. The line crosses z = 0 at most
once, so tau is exact, and each face is integrated with a composite Gauss-Legendre rule whose
panels meet at z = 0, where the integrand has a kink. The reference changes by less than 1e-11
when the panels are doubled.
"""

import importlib
import math
import os
import sys

if "opensn_console" not in globals():
    from mpi4py import MPI

    rank = MPI.COMM_WORLD.rank
    size = MPI.COMM_WORLD.size
    sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), "../../../../../")))
    from pyopensn.fieldfunc import FieldFunctionInterpolationVolume
    from pyopensn.logvol import RPPLogicalVolume
    from pyopensn.mesh import OrthogonalMeshGenerator
    from pyopensn.solver import UncollidedProblem, UncollidedSolver
    from pyopensn.source import PointSource
    from pyopensn.xs import MultiGroupXS

sys.path.append(os.path.dirname(__file__))
uncollided_utils = importlib.import_module("uncollided_unstructured_utils")
composite_rule = uncollided_utils.composite_rule
remove_file = uncollided_utils.remove_file
volume_integral = uncollided_utils.volume_integral

SIGMA_GROUND, SIGMA_ABOVE = 0.15, 1.0e-4
SOURCE = (1030.0, 1010.0, 100.0)
BOUNDS = [(0.0, 2000.0), (0.0, 2000.0), (-100.0, 2000.0)]
XY_NODES = [100.0 * i for i in range(21)]
Z_NODES = [100.0 * z for z in (-1.0, -0.6, -0.35, -0.2, -0.11, -0.06, -0.03, -0.01, 0.0, 0.4,
                               0.8, 1.4, 2.0, 3.0, 4.0, 6.0, 8.0, 11.0, 15.0, 20.0)]


def optical_depth(point):
    """Optical depth of the line from the source to a point on the box boundary."""
    length = math.dist(SOURCE, point)
    below = 0.0
    if point[2] < 0.0:
        below = length * (0.0 - point[2]) / (SOURCE[2] - point[2])
    return SIGMA_GROUND * below + SIGMA_ABOVE * (length - below)


def exact_escape(num_panels=24):
    total = 0.0
    for axis in range(3):
        for side in range(2):
            h = abs(BOUNDS[axis][side] - SOURCE[axis])
            a, b = [t for t in range(3) if t != axis]
            rules = []
            for t in (a, b):
                lo, hi = BOUNDS[t][0] - SOURCE[t], BOUNDS[t][1] - SOURCE[t]
                if t == 2 and lo < -SOURCE[2] < hi:
                    rules.append(composite_rule(lo, -SOURCE[2], num_panels)
                                 + composite_rule(-SOURCE[2], hi, num_panels))
                else:
                    rules.append(composite_rule(lo, hi, num_panels))
            for x, wx in rules[0]:
                for y, wy in rules[1]:
                    point = list(SOURCE)
                    point[axis] = BOUNDS[axis][side]
                    point[a] += x
                    point[b] += y
                    r = math.dist(SOURCE, point)
                    total += wx * wy * h / r**3 * math.exp(-optical_depth(point))
    return total / (4.0 * math.pi)


def outflow(exact_near_source_leakage):
    grid = OrthogonalMeshGenerator(node_sets=[XY_NODES, XY_NODES, Z_NODES]).Execute()
    grid.SetUniformBlockID(1)
    ground = RPPLogicalVolume(infx=True, infy=True, zmin=-200.0, zmax=0.0)
    above = RPPLogicalVolume(infx=True, infy=True, zmin=0.0, zmax=2100.0)
    grid.SetBlockIDFromLogicalVolume(ground, 0, True)
    xs_ground, xs_above = MultiGroupXS(), MultiGroupXS()
    xs_ground.CreateSimpleOneGroup(sigma_t=SIGMA_GROUND, c=0.0)
    xs_above.CreateSimpleOneGroup(sigma_t=SIGMA_ABOVE, c=0.0)

    file_name = f"uncollided_orthogonal_layered_rays_{int(exact_near_source_leakage)}.h5"
    remove_file(file_name)
    problem = UncollidedProblem(
        mesh=grid,
        num_groups=1,
        groupsets=[{"groups_from_to": [0, 0]}],
        xs_map=[{"block_ids": [0], "xs": xs_ground}, {"block_ids": [1], "xs": xs_above}],
        point_sources=[PointSource(location=list(SOURCE), strength=[1.0])],
        near_source=[RPPLogicalVolume(infx=True, infy=True, infz=True)],
        exact_near_source_leakage=exact_near_source_leakage,
        scattering_order=0,
    )
    solver = UncollidedSolver(problem=problem, file_name=file_name, progress_interval=0)
    solver.Initialize()
    solver.Execute()
    scalar_flux = problem.GetScalarFluxFieldFunction()[0]
    absorption = (SIGMA_GROUND * volume_integral(scalar_flux, ground,
                                                 FieldFunctionInterpolationVolume)
                  + SIGMA_ABOVE * volume_integral(scalar_flux, above,
                                                  FieldFunctionInterpolationVolume))
    remove_file(file_name)
    return 1.0 - absorption


if __name__ == "__main__":
    if size != 1:
        sys.exit(f"Expected one process, got {size}.")

    reference = exact_escape()
    error = abs(outflow(True) - reference) / reference
    if rank == 0:
        print(f"UncollidedOrthogonalLayeredOutflowError={error:.8e}")
