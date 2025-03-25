# Setup mesh
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

# Setting left, right, top and bottom boundaries
# left = 0
# right = 1
# bottom = 2
# top = 3
function dot_product(v1, v2)
  result = v1[1] * v2[1] + v1[2] * v2[2] + v1[3] * v2[3]
  return result
end

function bnd_id(pt, normal, cur_bid)
  epsilon = 1.0e-6
  n = { normal.x, normal.y, normal.z }
  if dot_product(n, { -1.0, 0.0, 0.0 }) > (1.0 - epsilon) then
    return 0
  end
  if dot_product(n, { 1.0, 0.0, 0.0 }) > (1.0 - epsilon) then
    return 1
  end
  if dot_product(n, { 0.0, -1.0, 0.0 }) > (1.0 - epsilon) then
    return 2
  end
  if dot_product(n, { 0.0, 1.0, 0.0 }) > (1.0 - epsilon) then
    return 3
  end

  return cur_bid
end

mesh.SetBoundaryIDFromFunction(grid, "bnd_id")

mesh.ExportToPVTU(grid, "new_bnd_ids")
