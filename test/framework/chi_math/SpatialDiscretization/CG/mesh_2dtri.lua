--############################################### Setup mesh
meshgen1 = chi_mesh.FromFileMeshGenerator.Create
({
  filename="../../../../../resources/TestMeshes/TriangleMesh2x2.obj"
})
chi_mesh.MeshGenerator.Execute(meshgen1)