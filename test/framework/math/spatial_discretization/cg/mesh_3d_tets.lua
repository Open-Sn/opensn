--############################################### Setup mesh
meshgen1 = mesh.FromFileMeshGenerator.Create({
  filename = "../../../../../resources/TestMeshes/GMSH_AllTets.vtu",
})
mesh.MeshGenerator.Execute(meshgen1)
