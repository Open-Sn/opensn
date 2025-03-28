--[[ @doc
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


--]]
-- Setup the mesh
meshgen = mesh.MeshGenerator.Create({
  inputs = {
    mesh.FromFileMeshGenerator.Create({
      filename = "../../../test/lua/modules/linear_boltzmann_solvers/transport_keigen/c5g7/mesh/2D_c5g7_coarse.msh",
    }),
  },
  partitioner = mesh.PETScGraphPartitioner.Create({ type = "parmetis" }),
})
grid = meshgen:Execute()

mesh.ExportToPVTU(grid, "c5g7_mesh_only")
