#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
A mesh that rank 0 rejects is rejected on every rank (2 ranks).

DistributedMeshGenerator reads and partitions the mesh on rank 0 only, while the other ranks wait
for their parts. The mesh of two tetrahedra on the same side of their shared face is folded, and
rank 0 rejects it while it builds the connectivity. Every rank must raise the error, so that no rank
is left waiting for mesh data. A valid mesh of two tetrahedra is then distributed, which completes
only if all the ranks are still in step.

OrthogonalMeshGenerator with distributed_generation=True also partitions on rank 0 only. A KBA cut
at x = 3 beyond the 2 x 2 mesh on [0, 2]^2 leaves the second slab without cells, so rank 0 rejects
the partition. Every rank must raise the error, instead of waiting for rank 0 in the exchange of
cell data. A valid distributed orthogonal mesh is then generated.
"""

import os
import sys

if "opensn_console" not in globals():
    from mpi4py import MPI

    rank = MPI.COMM_WORLD.rank
    global_sum = MPI.COMM_WORLD.allreduce
    sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), "../../../../")))
    from pyopensn.mesh import (
        DistributedMeshGenerator,
        FromFileMeshGenerator,
        KBAGraphPartitioner,
        OrthogonalMeshGenerator,
    )
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

    nodes = [0.0, 1.0, 2.0]
    try:
        OrthogonalMeshGenerator(
            node_sets=[nodes, nodes],
            distributed_generation=True,
            partitioner=KBAGraphPartitioner(nx=2, ny=1, xcuts=[3.0]),
        ).Execute()
        rejected = 0
    except RuntimeError:
        rejected = 1
    num_ortho_rejected = global_sum(rejected)

    ortho_grid = OrthogonalMeshGenerator(
        node_sets=[nodes, nodes],
        distributed_generation=True,
        partitioner=KBAGraphPartitioner(nx=2, ny=1, xcuts=[1.0]),
    ).Execute()
    num_ortho_cells = ortho_grid.GetGlobalNumberOfCells()

    if rank == 0:
        print(f"Ranks that rejected the folded mesh = {num_rejected}")
        print(f"Cells of the valid mesh = {num_cells}")
        print(f"Ranks that rejected the empty KBA slab = {num_ortho_rejected}")
        print(f"Cells of the valid distributed orthogonal mesh = {num_ortho_cells}")
