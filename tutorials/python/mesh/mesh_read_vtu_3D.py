#!/usr/bin/env python3

import os
import sys

if "opensn_console" not in globals():
    from mpi4py import MPI
    size = MPI.COMM_WORLD.size
    rank = MPI.COMM_WORLD.rank
    sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), "../../../../../")))
    from pyopensn.mesh import FromFileMeshGenerator, PETScGraphPartitioner

if __name__ == "__main__":
    """
    # Reading a 3D .vtu file

    ## Read the Mesh
    We start by reading a 3D vtu file. The resulting mesh and partition is shown below:

    ![Mesh_Partition](images/read_3d_vtu_partition.png)
    """
    # Setup mesh
    meshgen = FromFileMeshGenerator(
        filename="../../../test/assets/mesh/GMSH_AllTets.vtu",
        partitioner=PETScGraphPartitioner(type='parmetis'),
    )
    grid = meshgen.Execute()

    # Set block IDs
    grid.SetUniformBlockID(0)

    # Export
    grid.ExportToPVTU("Read_3D_mesh_only")
