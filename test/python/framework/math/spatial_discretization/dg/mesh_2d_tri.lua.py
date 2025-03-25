# Setup mesh
meshgen1 = mesh.FromFileMeshGenerator.Create({
  filename = "+/+/+/+/+/assets/mesh/TriangleMesh2x2.obj",
})
grid = meshgen1:Execute()
