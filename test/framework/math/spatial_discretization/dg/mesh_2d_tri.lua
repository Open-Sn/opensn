--############################################### Setup mesh
meshgen1 = mesh.FromFileMeshGenerator.Create({
  filename = "../../../../assets/mesh/TriangleMesh2x2.obj",
})
mesh.MeshGenerator.Execute(meshgen1)
