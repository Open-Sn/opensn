#!/usr/bin/env python3
"""Verify uncollided generation with optically thick cells at the vacuum boundary.

A point source sits near the center of a 20 cm cube of a pure absorber (sigma_t = 1/cm). The mesh
is refined around the source and coarsens to cells 3 cm (3 mean free paths) thick at the
boundary. In those cells the PWLD flux is negative on some outgoing faces, so the net outflow is
slightly negative and the removal slightly exceeds the production, while the balance still holds.
Generation used to fail on this with "physical uncollided removal exceeds the physical source
rate". It now completes and warns.

The removal, sigma_t times the volume integral of the scalar flux, is compared with one minus the
exact escape fraction of the box (1.38e-5). The PWLD outflow error at the thick boundary cells
is 2.0e-5 of the source.

Checks:
  UncollidedThickBoundaryNegativeOutflow: 1 if 1 - removal is negative, so the case is exercised.
  UncollidedThickBoundaryOutflowError: |(1 - removal) - exact escape fraction| (below 1e-4).
"""

import importlib
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
remove_file = uncollided_utils.remove_file
uncollided_escape_3d_box = uncollided_utils.uncollided_escape_3d_box
volume_integral = uncollided_utils.volume_integral


if __name__ == "__main__":
    if size != 1:
        sys.exit(f"Expected one process, got {size}.")

    nodes = [-6.0, -3.0, -1.0, 0.0, 1.0] + [2.0 + 0.4 * i for i in range(11)]
    nodes += [7.0, 8.0, 9.0, 11.0, 14.0]
    grid = OrthogonalMeshGenerator(node_sets=[nodes, nodes, nodes]).Execute()
    grid.SetUniformBlockID(0)

    sigma_t = 1.0
    xs = MultiGroupXS()
    xs.CreateSimpleOneGroup(sigma_t=sigma_t, c=0.0)

    source = (3.03, 4.11, 3.93)
    near_source = RPPLogicalVolume(
        xmin=source[0] - 0.6,
        xmax=source[0] + 0.6,
        ymin=source[1] - 0.6,
        ymax=source[1] + 0.6,
        zmin=source[2] - 0.6,
        zmax=source[2] + 0.6,
    )
    file_name = "uncollided_thick_boundary_outflow.h5"
    remove_file(file_name)
    problem = UncollidedProblem(
        mesh=grid,
        num_groups=1,
        groupsets=[{"groups_from_to": [0, 0]}],
        xs_map=[{"block_ids": [0], "xs": xs}],
        point_sources=[PointSource(location=list(source), strength=[1.0])],
        near_source=[near_source],
        scattering_order=0,
    )
    solver = UncollidedSolver(problem=problem, file_name=file_name)
    solver.Initialize()
    solver.Execute()
    remove_file(file_name)

    whole_domain = RPPLogicalVolume(
        xmin=-6.0, xmax=14.0, ymin=-6.0, ymax=14.0, zmin=-6.0, zmax=14.0
    )
    removal = sigma_t * volume_integral(
        problem.GetScalarFluxFieldFunction()[0], whole_domain, FieldFunctionInterpolationVolume
    )
    outflow = 1.0 - removal
    exact = uncollided_escape_3d_box(sigma_t, source, [(-6.0, 14.0)] * 3)

    if rank == 0:
        print(f"UncollidedThickBoundaryOutflow={outflow:.8e}")
        print(f"UncollidedThickBoundaryExactEscape={exact:.8e}")
        print(f"UncollidedThickBoundaryNegativeOutflow={int(outflow < 0.0)}")
        print(f"UncollidedThickBoundaryOutflowError={abs(outflow - exact):.8e}")
        sys.stdout.flush()
