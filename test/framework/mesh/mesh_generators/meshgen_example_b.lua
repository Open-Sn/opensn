meshgen1 = mesh.FromFileMeshGenerator.Create
({
  filename="triangle_mesh_2x2.obj"
})
mesh.MeshGenerator.Execute(meshgen1)

--MeshHandlerExportMeshToVTK("ZMeshTest")
