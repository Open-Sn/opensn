#!/usr/bin/env python3
"""Verify that uncollided generation rejects multiple MPI ranks."""

import os
import sys

if "opensn_console" not in globals():
    from mpi4py import MPI

    rank = MPI.COMM_WORLD.rank
    size = MPI.COMM_WORLD.size
    sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), "../../../../../")))
    from pyopensn.logvol import RPPLogicalVolume
    from pyopensn.mesh import FromFileMeshGenerator
    from pyopensn.solver import UncollidedProblem, UncollidedSolver
    from pyopensn.source import PointSource
    from pyopensn.xs import MultiGroupXS

from uncollided_unstructured_utils import mesh_path


if __name__ == "__main__":
    if size != 2:
        sys.exit(f"Expected two processes, got {size}.")

    grid = FromFileMeshGenerator(filename=mesh_path("triangle_mesh2x2_fine.obj")).Execute()
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
    rejected = False
    try:
        problem = UncollidedProblem(
            mesh=grid,
            num_groups=1,
            groupsets=[{"groups_from_to": [0, 0]}],
            xs_map=[{"block_ids": [0], "xs": xs}],
            point_sources=[PointSource(location=[0.0, 0.0, 0.0], strength=[1.0])],
            near_source=[whole_domain],
        )
        solver = UncollidedSolver(
            problem=problem,
            file_name="uncollided_parallel_rejected.h5",
        )
        solver.Initialize()
    except ValueError as error:
        rejected = "exactly one MPI rank" in str(error)

    if rank == 0:
        print(f"UncollidedParallelGenerationRejected={int(rejected)}")
        sys.stdout.flush()

    if not rejected:
        raise RuntimeError("Uncollided generation did not reject multiple MPI ranks")
