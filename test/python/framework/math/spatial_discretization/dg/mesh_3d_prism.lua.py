# Setup mesh
meshgen1 = mesh.ExtruderMeshGenerator.Create({
  inputs = {
    mesh.FromFileMeshGenerator.Create({
      filename = "+/+/+/+/+/assets/mesh/TriangleMesh2x2.obj",
    }),
  },
  layers = { { z = 0.4, n = 2 }, { z = 0.8, n = 2 }, { z = 1.2, n = 2 }, { z = 1.6, n = 2 } }, # layers
})
grid = meshgen1:Execute()
