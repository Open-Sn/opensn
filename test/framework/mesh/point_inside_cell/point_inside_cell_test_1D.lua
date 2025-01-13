nodes = { -1.0, -0.75, 0.0, 1.0, 2.0 }
meshgen = mesh.OrthogonalMeshGenerator.Create({
  node_sets = { nodes }
})
mesh.MeshGenerator.Execute(meshgen)

unit_tests.point_inside_cell_test_00()
