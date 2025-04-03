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
    # 2D Orthogonal Grid
    Here, we will use the in-house orthogonal mesh generator for a simple Cartesian grid.

    ## List of nodes
    We first create a list of nodes. The nodes will be spread from -L/2 to +L/2.

    ## Mesh and KBA partition
    We use the `OrthogonalMeshGenerator` and pass the list of nodes per dimension. Here, we pass
    2 times the same list of nodes to create a 2D geometry with square cells. Thus, we create a
    square domain, of side length L, centered on the origin (0,0).

    We also introduce partitioning. The following schemes are available
    - KBA partitioning for regular grids that can be cut into right parallelpipeds
    - Parmetis partitioning and Scotch partitioning, applicable with any grid type.

    Here, we first partition the 2D mesh into 2x2 subdomains using `KBAGraphPartitioner`. We place
    one cut along the x-axis at x=0 by filling the xcuts array. Likewise for ycuts (y=0).
    The cell assignment to a partition is done based
    on where the cell center is located with respect to the various xcuts, ycuts, and zcuts
    (a fuzzy logic is applied to avoid issues).

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
        node_sets=[nodes, nodes],
        partitioner=KBAGraphPartitioner(
            nx=2,
            ny=2,
            xcuts=[0.0],
            ycuts=[0.0],
        ),
    )
    grid = meshgen.Execute()

    # Set block IDs
    grid.SetUniformBlockID(0)

    grid.ExportToPVTU("ortho_2D_KBA")

    meshgen2 = OrthogonalMeshGenerator(
        node_sets=[nodes, nodes],
        partitioner=PETScGraphPartitioner(type='parmetis'),
    )
    grid2 = meshgen2.Execute()

    # Set block IDs
    grid2.SetUniformBlockID(0)

    grid2.ExportToPVTU("ortho_2D_Parmetis")
