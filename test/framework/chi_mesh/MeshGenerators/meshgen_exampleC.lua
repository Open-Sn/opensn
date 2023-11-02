meshgen1 = chi_mesh.ExtruderMeshGenerator.Create
({
  inputs =
  {
    chi_mesh.FromFileMeshGenerator.Create
    ({
      filename="TriangleMesh2x2.obj"
    }),
  },
  layers = {{z=1.1, n=2}, {z=2.1, n=3}}
})
chi_mesh.MeshGenerator.Execute(meshgen1)

--chiMeshHandlerExportMeshToVTK("ZMeshTest")
