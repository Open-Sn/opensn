--[[ @doc
# Reading a 3D .vtu file

## Read the Mesh
We start by reading a 3D vtu file. The resulting mesh and partition is shown below:

![Mesh_Partition](images/read_3d_vtu_partition.png)

--]]
-- Setup the mesh
meshgen = mesh.MeshGenerator.Create({
  inputs = {
    mesh.FromFileMeshGenerator.Create({
      filename = "../../test/assets/mesh/GMSH_AllTets.vtu",
    }),
  },
  partitioner = mesh.PETScGraphPartitioner.Create({ type = "parmetis" }),
})
meshgen:Execute()

mesh.ExportToPVTU("Read_3D_mesh_only")
