x_nodes = { -1.0, -0.75, -0.5, -0.25, 0.0, 0.5, 1.0 }
y_nodes = { -1.0, -0.75, -0.5, -0.25, 0.0, 0.25, 0.5, 0.75, 1.0 }
fine = mesh.OrthogonalMeshGenerator.Create({
  node_sets = { x_nodes, y_nodes }
})
mesh.MeshGenerator.Execute(fine)

nodes = { -1.0, -0.5, 0.0, 0.5, 1.0 }
coarse = mesh.OrthogonalMeshGenerator.Create({
  node_sets = { nodes, nodes }
})
mesh.MeshGenerator.Execute(coarse)

mesh.MeshMapping.Build()
