-- Setup mesh
meshgen1 = mesh.FromFileMeshGenerator.Create({
  filename = "../../../../assets/mesh/TriangleMesh2x2.obj",
})
meshgen1:Execute()
