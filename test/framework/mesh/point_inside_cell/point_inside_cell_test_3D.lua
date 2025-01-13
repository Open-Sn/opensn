x_nodes = { -1.0, 1.0 }
y_nodes = { 0.0, 0.5, 1.0 }
z_nodes = { -1.0, 0.0, 1.0 }
meshgen = mesh.OrthogonalMeshGenerator.Create({
  node_sets = { x_nodes, y_nodes, z_nodes }
})
mesh.MeshGenerator.Execute(meshgen)

unit_tests.point_inside_cell_test_00()
