-- Setup mesh
nodes = {}
N = 10
L = 5
xmin = -L / 2
dx = L / N
for i = 1, (N + 1) do
  k = i - 1
  nodes[i] = xmin + k * dx
end

meshgen1 = mesh.OrthogonalMeshGenerator.Create({ node_sets = { nodes, nodes, nodes } })
grid = meshgen1:Execute()

--Sets a middle square to material 1
function mat_id(pt, cur_id)
  if math.abs(pt.x) < L / 10 and math.abs(pt.y) < L / 10 then
    return 1
  end
  return cur_id
end

mesh.SetBlockIDFromFunction(grid, "mat_id")

mesh.ExportToPVTU(grid, "new_mat_ids")
