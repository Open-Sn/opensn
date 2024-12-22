fine = { -1.0, -0.75, -0.5, -0.25, 0.0, 0.5, 1.0 }
meshgen1 = mesh.OrthogonalMeshGenerator.Create({
  node_sets = { nodes }
})
mesh.MeshGenerator.Execute(fine)

nodes = { -1.0, -0.5, 0.0, 0.5, 1.0 }
coarse = mesh.OrthogonalMeshGenerator.Create({
  node_sets = { nodes }
})
mesh.MeshGenerator.Execute(coarse)

mesh.MeshMapping.Build()
