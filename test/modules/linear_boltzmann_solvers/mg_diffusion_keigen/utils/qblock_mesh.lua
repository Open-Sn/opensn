--############################################### Setup mesh
nodes = {}
N = 40
L = 14.0
xmin = 0.0
dx = L / N
for i = 1, (N + 1) do
  k = i - 1
  nodes[i] = xmin + k * dx
end

meshgen1 = mesh.OrthogonalMeshGenerator.Create({ node_sets = { nodes, nodes } })
mesh.MeshGenerator.Execute(meshgen1)

mesh.SetUniformMaterialID(0)

vol1 =
  logvol.RPPLogicalVolume.Create({ xmin = -1000.0, xmax = 10.0, ymin = -1000.0, ymax = 10.0, infz = true })
mesh.SetMaterialIDFromLogicalVolume(vol1, 1)
