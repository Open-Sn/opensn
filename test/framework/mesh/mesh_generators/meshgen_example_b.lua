meshgen1 = chi_mesh.FromFileMeshGenerator.Create
({
  filename="triangle_mesh_2x2.obj"
})
chi_mesh.MeshGenerator.Execute(meshgen1)

--chiMeshHandlerExportMeshToVTK("ZMeshTest")
