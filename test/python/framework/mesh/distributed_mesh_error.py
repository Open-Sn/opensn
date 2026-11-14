#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
A mesh that rank 0 rejects is rejected on every rank (2 ranks).

DistributedMeshGenerator reads and partitions the mesh on rank 0 only, while the other ranks wait
for their parts. The mesh of two tetrahedra on the same side of their shared face is folded, and
rank 0 rejects it while it builds the connectivity. Every rank must raise the error, so that no rank
is left waiting for mesh data. A valid mesh of two tetrahedra is then distributed, which completes
only if all the ranks are still in step.
"""

import os
import sys

if "opensn_console" not in globals():
    from mpi4py import MPI

    rank = MPI.COMM_WORLD.rank
    global_sum = MPI.COMM_WORLD.allreduce
    sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), "../../../../")))
    from pyopensn.mesh import DistributedMeshGenerator, FromFileMeshGenerator
else:
    global_sum = MPIAllReduce

if __name__ == "__main__":
    mesh_dir = "../../../assets/mesh"

    try:
        DistributedMeshGenerator(
            inputs=[FromFileMeshGenerator(filename=os.path.join(mesh_dir, "two_tets_folded.msh"))]
        ).Execute()
        rejected = 0
    except RuntimeError:
        rejected = 1
    num_rejected = global_sum(rejected)

    grid = DistributedMeshGenerator(
        inputs=[FromFileMeshGenerator(filename=os.path.join(mesh_dir, "two_tets_valid.msh"))]
    ).Execute()
    num_cells = grid.GetGlobalNumberOfCells()

    if rank == 0:
        print(f"Ranks that rejected the folded mesh = {num_rejected}")
        print(f"Cells of the valid mesh = {num_cells}")
