nodes = {-1.0,-0.75,-0.5,-0.25,0.0,0.25,0.5,0.75,1.0}
meshgen1 = mesh.OrthogonalMeshGenerator.Create
({
  node_sets = {nodes,nodes},
})
mesh.MeshGenerator.Execute(meshgen1)
