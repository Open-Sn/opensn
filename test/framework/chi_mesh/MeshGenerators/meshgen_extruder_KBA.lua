meshgen1 = chi_mesh.ExtruderMeshGenerator.Create
({
  inputs =
  {
    chi_mesh.FromFileMeshGenerator.Create
    ({
      filename="TriangleMesh2x2.obj"
    }),
  },
  layers = {{z=1.1, n=2}, {z=2.1, n=3}},
  partitioner = chi.KBAGraphPartitioner.Create
  ({
    nx = 2, ny=2, nz=2,
    xcuts = {0.0}, ycuts = {0.0}, zcuts = {1.1}
  })
})
chi_mesh.MeshGenerator.Execute(meshgen1)

--chiMeshHandlerExportMeshToVTK("ZMeshTest")
