nodes = {}
N = 40
L = 5
xmin = -L / 2
dx = L / N
for i = 1, (N + 1) do
  k = i - 1
  nodes[i] = xmin + k * dx
end
meshgen1 = mesh.OrthogonalMeshGenerator.Create({
  node_sets = { nodes, nodes, nodes },
})
grid = meshgen1:Execute()

lv1 = logvol.RCCLogicalVolume.Create({ r = 1.3, x0 = L / 2, y0 = L / 2, z0 = -1.0, vz = 2.0 })
grid:SetBlockIDFromLogicalVolume(lv1, 1, true)

lv2 = logvol.RCCLogicalVolume.Create({
  r = 1.3,
  x0 = -0.8,
  y0 = -0.8,
  z0 = -1.5,
  vx = 1.0,
  vy = 1.0,
  vz = 3.0,
})
grid:SetBlockIDFromLogicalVolume(lv2, 2, true)

mesh.ExportToPVTU(grid, "lv_rcc_test1")
