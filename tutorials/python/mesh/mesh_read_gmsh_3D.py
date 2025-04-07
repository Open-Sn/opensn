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
    # Reading a 3D .msh file

    Reading a 3D .msh file generated with the [gmsh Mesh Generator](https://gmsh.info/).

    ## Read the Mesh
    The resulting mesh and material layout is shown below:

    ![Mesh_Material](images/c5g7_coarse_material.png)

    When using the Parmetis partitioner, we obtain:

    ![Mesh_Partition](images/c5g7_coarse_partition.png)
    """
    # Setup mesh
    meshgen = FromFileMeshGenerator(
        filename="./two_spheres_small.msh",
        partitioner=PETScGraphPartitioner(type='parmetis'),
    )
    grid = meshgen.Execute()

    # Export
    grid.ExportToPVTU("Read_3D_mesh_only")
