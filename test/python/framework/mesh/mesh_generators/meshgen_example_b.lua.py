meshgen1 = mesh.FromFileMeshGenerator.Create({
  filename = "triangle_mesh_2x2.obj",
})
grid = meshgen1:Execute()
