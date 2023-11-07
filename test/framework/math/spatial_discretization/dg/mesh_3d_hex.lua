--############################################### Setup mesh
meshgen1 = chi_mesh.FromFileMeshGenerator.Create
({
  filename="../../../../../resources/TestMeshes/GMSH_AllHexes.vtu"
})
chi_mesh.MeshGenerator.Execute(meshgen1)