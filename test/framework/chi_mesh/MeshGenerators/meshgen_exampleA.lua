nodes = {-1.0,-0.75,-0.5,-0.25,0.0,0.25,0.5,0.75,1.0}
meshgen1 = chi_mesh.OrthogonalMeshGenerator.Create
({
  node_sets = {nodes,nodes},
})
chi_mesh.MeshGenerator.Execute(meshgen1)

--chiMeshHandlerExportMeshToVTK("ZMeshTest")
