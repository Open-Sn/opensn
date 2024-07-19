--############################################### Setup mesh
meshgen1 = mesh.FromFileMeshGenerator.Create({
  filename = "../../../../assets/mesh/GMSH_AllHexes.vtu",
})
mesh.MeshGenerator.Execute(meshgen1)
