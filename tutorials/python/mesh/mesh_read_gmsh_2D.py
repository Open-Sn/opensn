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
    # Reading a 2D .msh file

    Reading a 2D .msh file with material IDs and boundary IDs.

    We have created an unstructured mesh with the [gmsh Mesh Generator](https://gmsh.info/).

    ## Read the Mesh
    We use the ```FromFileMeshGenerator``` and pass the path to the msh file.
    We also partition the 2D mesh into 2x2 subdomains using `Parmetis`.
    Finally, we export the mesh to a VTU file.

    The resulting mesh and material layout is shown below:

    ![Mesh_Material](images/c5g7_coarse_material.png)

    When using the Parmetis partitioner, we obtain:

    ![Mesh_Partition](images/c5g7_coarse_partition.png)
    """
    # Setup mesh
    filepath = "../../../test/python/modules/linear_boltzmann_solvers/transport_keigen/c5g7/"\
               + "mesh/2D_c5g7_coarse.msh"
    meshgen = FromFileMeshGenerator(
        filename=filepath,
        partitioner=PETScGraphPartitioner(type='parmetis'),
    )
    grid = meshgen.Execute()

    # Export
    grid.ExportToPVTU("c5g7_mesh_only")
