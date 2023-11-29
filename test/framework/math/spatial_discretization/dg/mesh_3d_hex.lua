--############################################### Setup mesh
meshgen1 = mesh.FromFileMeshGenerator.Create
({
  filename="../../../../../resources/TestMeshes/GMSH_AllHexes.vtu"
})
mesh.MeshGenerator.Execute(meshgen1)
