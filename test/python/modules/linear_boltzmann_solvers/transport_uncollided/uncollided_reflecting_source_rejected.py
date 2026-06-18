#!/usr/bin/env python3
"""Verify that reflecting uncollided setups reject face- and vertex-aligned sources."""

import os
import sys

if "opensn_console" not in globals():
    from mpi4py import MPI

    rank = MPI.COMM_WORLD.rank
    size = MPI.COMM_WORLD.size
    sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), "../../../../../")))
    from pyopensn.logvol import RPPLogicalVolume
    from pyopensn.mesh import OrthogonalMeshGenerator
    from pyopensn.solver import UncollidedProblem
    from pyopensn.source import PointSource
    from pyopensn.xs import MultiGroupXS


if __name__ == "__main__":
    if size != 1:
        sys.exit(f"Expected one process, got {size}.")

    grid = OrthogonalMeshGenerator(node_sets=[[-1.0, 0.0, 1.0], [-1.0, 0.0, 1.0]]).Execute()
    grid.SetUniformBlockID(0)
    grid.SetOrthogonalBoundaries()

    xs = MultiGroupXS()
    xs.CreateSimpleOneGroup(sigma_t=1.0, c=0.0)

    whole_domain = RPPLogicalVolume(
        xmin=-1.01,
        xmax=1.01,
        ymin=-1.01,
        ymax=1.01,
        infz=True,
    )

    rejected = False
    try:
        UncollidedProblem(
            mesh=grid,
            num_groups=1,
            groupsets=[{"groups_from_to": [0, 0]}],
            xs_map=[{"block_ids": [0], "xs": xs}],
            boundary_conditions=[{"name": "xmin", "type": "reflecting"}],
            point_sources=[PointSource(location=[0.0, 0.25, 0.0], strength=[1.0])],
            near_source=[whole_domain],
        )
    except Exception as error:
        rejected = "not currently supported for uncollided generation" in str(error)

    if rank == 0:
        print(f"UncollidedReflectingFaceOrVertexSourceRejected={int(rejected)}")
        sys.stdout.flush()

    if not rejected:
        raise RuntimeError(
            "Reflecting uncollided generation accepted a point source on a face, edge, or vertex."
        )
