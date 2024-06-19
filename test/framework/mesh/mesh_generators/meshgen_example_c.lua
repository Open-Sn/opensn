meshgen1 = mesh.ExtruderMeshGenerator.Create({
  inputs = {
    mesh.FromFileMeshGenerator.Create({
      filename = "triangle_mesh_2x2.obj",
    }),
  },
  layers = { { z = 1.1, n = 2 }, { z = 2.1, n = 3 } },
})
mesh.MeshGenerator.Execute(meshgen1)
