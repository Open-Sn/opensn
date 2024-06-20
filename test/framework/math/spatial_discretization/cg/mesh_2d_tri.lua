--############################################### Setup mesh
meshgen1 = mesh.FromFileMeshGenerator.Create({
  filename = "../../../../../resources/TestMeshes/TriangleMesh2x2.obj",
})
mesh.MeshGenerator.Execute(meshgen1)
