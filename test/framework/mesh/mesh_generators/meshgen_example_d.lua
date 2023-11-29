meshgen1 = mesh.MeshGenerator.Create
({
  inputs =
  {
    mesh.FromFileMeshGenerator.Create
    ({
      filename="triangle_mesh_2x2.obj"
    }),
    mesh.ExtruderMeshGenerator.Create
    ({
      layers = {{z=1.1, n=2}, {z=2.1, n=3}}
    })
  }
})
mesh.MeshGenerator.Execute(meshgen1)

--chiMeshHandlerExportMeshToVTK("ZMeshTest")
