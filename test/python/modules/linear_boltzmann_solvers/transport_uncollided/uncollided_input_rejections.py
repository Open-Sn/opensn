#!/usr/bin/env python3
"""Verify common uncollided input-validation failures."""

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

    xs = MultiGroupXS()
    xs.CreateSimpleOneGroup(sigma_t=1.0, c=0.0)

    whole_domain = RPPLogicalVolume(
        xmin=-1.01,
        xmax=1.01,
        ymin=-1.01,
        ymax=1.01,
        infz=True,
    )

    missing_point_source_rejected = False
    try:
        UncollidedProblem(
            mesh=grid,
            num_groups=1,
            groupsets=[{"groups_from_to": [0, 0]}],
            xs_map=[{"block_ids": [0], "xs": xs}],
            near_source=[whole_domain],
        )
    except Exception as error:
        missing_point_source_rejected = "at least one point source is required" in str(error)

    near_source_mismatch_rejected = False
    try:
        UncollidedProblem(
            mesh=grid,
            num_groups=1,
            groupsets=[{"groups_from_to": [0, 0]}],
            xs_map=[{"block_ids": [0], "xs": xs}],
            point_sources=[PointSource(location=[0.25, -0.2, 0.0], strength=[1.0])],
            near_source=[],
        )
    except Exception as error:
        near_source_mismatch_rejected = (
            "number of near-source logical volumes must match the number of point sources"
            in str(error)
        )

    face_or_vertex_source_rejected = False
    try:
        UncollidedProblem(
            mesh=grid,
            num_groups=1,
            groupsets=[{"groups_from_to": [0, 0]}],
            xs_map=[{"block_ids": [0], "xs": xs}],
            point_sources=[PointSource(location=[0.0, 0.25, 0.0], strength=[1.0])],
            near_source=[whole_domain],
        )
    except Exception as error:
        face_or_vertex_source_rejected = (
            "not currently supported for uncollided generation" in str(error)
        )

    if rank == 0:
        print(f"UncollidedMissingPointSourceRejected={int(missing_point_source_rejected)}")
        print(f"UncollidedNearSourceMismatchRejected={int(near_source_mismatch_rejected)}")
        print(f"UncollidedFaceOrVertexSourceRejected={int(face_or_vertex_source_rejected)}")
        sys.stdout.flush()

    if not missing_point_source_rejected:
        raise RuntimeError("Uncollided generation accepted an input with no point source.")
    if not near_source_mismatch_rejected:
        raise RuntimeError("Uncollided generation accepted a near-source count mismatch.")
    if not face_or_vertex_source_rejected:
        raise RuntimeError(
            "Uncollided generation accepted a point source on a face, edge, or vertex."
        )
