-- Setup mesh

meshgen1 = mesh.MeshGenerator.Create({
  inputs = {
    mesh.FromFileMeshGenerator.Create({
      filename = "mesh/2D_c5g7_coarse.msh",
    }),
  },
})
grid = meshgen1:Execute()
