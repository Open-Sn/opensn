-- Setup mesh
meshgen1 = mesh.FromFileMeshGenerator.Create({
  filename = "../../../../assets/mesh/GMSH_AllTets.vtu",
})
meshgen1:Execute()
