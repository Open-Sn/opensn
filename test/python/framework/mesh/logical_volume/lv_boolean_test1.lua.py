--[[
test for boolean operations on logical volumes:
  lv1 = sphere
  lv2 = right circular cylinder (rcc)
  lv3 = in rcc but not in sphere
--]]

-- set up orthogonal 3D geometry
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

-- assign block ID 10 to all cells
vol0 = logvol.RPPLogicalVolume.Create({ infx = true, infy = true, infz = true })
grid:SetBlockIDFromLogicalVolume(vol0, 10, true)

-- create logical volume lv1 as an analytical sphere
lv1 = logvol.SphereLogicalVolume.Create({ r = 1.3, x = 1.0, y = -1.0, z = 2.0 })

-- create logical volume lv2 as an analytical rcc
lv2 = logvol.RCCLogicalVolume.Create({
  r = 1.3,
  x0 = -0.8,
  y0 = -0.8,
  z0 = -1.5,
  vx = 1.0,
  vy = 1.0,
  vz = 3.0,
})

-- create logical volume lv3 as boolean: true if cell is in lv2 and false if in lv1
lv3 = logvol.BooleanLogicalVolume.Create({
  parts = { { op = true, lv = lv2 }, { op = false, lv = lv1 } },
})

-- assign block ID 1 to all cells in lv3 which is the part of lv2 that is not in lv1
grid:SetBlockIDFromLogicalVolume(lv3, 1, true)

-- assign block ID 5 to all cells in lv1
grid:SetBlockIDFromLogicalVolume(lv1, 5, true)

-- export to vtk
mesh.ExportToPVTU(grid, "lv_boolean_test1")
