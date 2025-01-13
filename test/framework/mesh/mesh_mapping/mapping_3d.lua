xz_nodes = { -1.0, -0.75, -0.5, -0.25, 0.0, 0.5, 1.0 }
y_nodes = { -1.0, -0.75, -0.5, -0.25, 0.0, 0.25, 0.5, 0.75, 1.0 }
fine = mesh.OrthogonalMeshGenerator.Create({
  node_sets = { xz_nodes, y_nodes, xz_nodes }
  -- node_sets = { xz_nodes, yz_nodes, xz_nodes }
})
mesh.MeshGenerator.Execute(fine)

xy_nodes = { -1.0, 0.0,  1.0 }
z_nodes = { -1.0, -0.5, 0.0, 0.5, 1.0 }
coarse = mesh.OrthogonalMeshGenerator.Create({
  node_sets = { xy_nodes, xy_nodes, z_nodes }
})
mesh.MeshGenerator.Execute(coarse)

mesh.MeshMapping.Build()
