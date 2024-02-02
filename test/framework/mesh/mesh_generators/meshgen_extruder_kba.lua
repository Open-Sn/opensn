meshgen1 = mesh.ExtruderMeshGenerator.Create
({
  inputs =
  {
    mesh.FromFileMeshGenerator.Create
    ({
      filename="triangle_mesh_2x2.obj"
    }),
  },
  layers = {{z=1.1, n=2}, {z=2.1, n=3}},
  partitioner = KBAGraphPartitioner.Create
  ({
    nx = 2, ny=2, nz=2,
    xcuts = {0.0}, ycuts = {0.0}, zcuts = {1.1}
  })
})
mesh.MeshGenerator.Execute(meshgen1)

--MeshHandlerExportMeshToVTK("ZMeshTest")
