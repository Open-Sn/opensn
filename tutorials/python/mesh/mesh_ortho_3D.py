#!/usr/bin/env python3

import os
import sys

if "opensn_console" not in globals():
    from mpi4py import MPI
    size = MPI.COMM_WORLD.size
    rank = MPI.COMM_WORLD.rank
    sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), "../../../../../")))
    from pyopensn.mesh import OrthogonalMeshGenerator, KBAGraphPartitioner, PETScGraphPartitioner

if __name__ == "__main__":
    """
    # 3D Orthogonal Grid

    This is very similar to the previous `2D Orthogonal Grid` tutorial. We use the two
    partitioners again (KBA and Parmetis).

    ## Mesh and KBA partition
    A simple orthogonal 3D mesh with KBA partitioner. We will run with 8 processes, so we will
    partition the domain in 2x2x2 parts with cuts placed exactly at x=0, y=0, and z=0.

    ## Material IDs
    When using the in-house `OrthogonalMeshGenerator`, no material IDs are assigned. The user needs
    to assign material IDs to all cells. Here, we have a homogeneous domain, so we assign a material
    ID with value 0 for each cell in the spatial domain.

    ## Export the mesh
    We export to vtu format. The resulting mesh partition is shown below

    ## Mesh (again) and Parmetis partition
    Now, we use the Parmetis partitioner.

    """
    # Setup mesh
    length = 2.
    n_cells = 10
    dx = length / n_cells
    nodes = [i * dx for i in range(n_cells + 1)]

    # Setup mesh
    meshgen = OrthogonalMeshGenerator(
        node_sets=[nodes, nodes, nodes],
        partitioner=KBAGraphPartitioner(
            nx=2,
            ny=2,
            nz=2,
            xcuts=[0.0],
            ycuts=[0.0],
            zcuts=[0.0],
        ),
    )
    grid = meshgen.Execute()

    # Set block IDs
    grid.SetUniformBlockID(0)

    grid.ExportToPVTU("ortho_3D_KBA")

    meshgen2 = OrthogonalMeshGenerator(
        node_sets=[nodes, nodes],
        partitioner=PETScGraphPartitioner(type='parmetis'),
    )
    grid2 = meshgen2.Execute()

    # Set block IDs
    grid2.SetUniformBlockID(0)

    grid2.ExportToPVTU("ortho_3D_Parmetis")
