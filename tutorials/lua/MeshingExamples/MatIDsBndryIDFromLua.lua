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

mesh.ExportToPVTU("ZMeshPhase1")

grid:SetUniformBlockID(0)

--Sets a middle square to material 1
function MatIDFunction1(pt, cur_id)
  if math.abs(pt.x) < L / 10 and math.abs(pt.y) < L / 10 then
    return 1
  end

  return cur_id
end

mesh.SetBlockIDFromFunction("MatIDFunction1")

mesh.ExportToPVTU("ZMeshPhase2")

--Setting left, right, top and bottom boundaries
-- left = 0
-- right = 1
-- bottom = 2
-- top = 3
function dot_product(v1, v2)
  result = v1[1] * v2[1] + v1[2] * v2[2] + v1[3] * v2[3]
  return result
end

function BndryIDFunction1(pt, normal, cur_bid)
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

mesh.SetBoundaryIDFromFunction("BndryIDFunction1")

mesh.ExportToPVTU("ZMeshPhase3")
