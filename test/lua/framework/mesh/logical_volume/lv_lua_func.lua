-- test for lua function used as a logical volume
-- set up orthogonal 3D geometry
nodes = {}
N = 50
L = 5.0
xmin = -L / 2
dx = L / N
for i = 1, (N + 1) do
  k = i - 1
  nodes[i] = xmin + k * dx
end

meshgen = mesh.OrthogonalMeshGenerator.Create({
  node_sets = { nodes, nodes, nodes },
})
grid = meshgen:Execute()

-- assign mat ID 10 to whole domain
vol0 = logvol.RPPLogicalVolume.Create({ infx = true, infy = true, infz = true })
grid:SetBlockIDFromLogicalVolume(vol0, 10, true)

--Sets lua function describing a sphere (material 11)
function MatIDFunction1(pt, cur_id)
  if pt.x * pt.x + pt.y * pt.y + pt.z * pt.z < 1.0 then
    return 11
  end
  return cur_id
end

-- assign block ID 11 to lv using lua function
mesh.SetBlockIDFromFunction(grid, "MatIDFunction1")

-- export to vtk
mesh.ExportToPVTU(grid, "lv_lua_func_out")
