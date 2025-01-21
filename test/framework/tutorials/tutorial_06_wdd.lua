-- Setup mesh
if nmesh == nil then
  nmesh = 10
end

nodes = {}
N = nmesh
L = 2.0
xmin = -L / 2
dx = L / N
for i = 1, (N + 1) do
  k = i - 1
  nodes[i] = xmin + k * dx
end

meshgen1 = mesh.OrthogonalMeshGenerator.Create({ node_sets = { nodes, nodes } })
grid = meshgen1:Execute()

-- Set block IDs
grid:SetUniformBlockID(0)

unit_sim_tests.SimTest06_WDD({ mesh = grid })
MPIBarrier()
if location_id == 0 then
  os.execute("rm SimTest_06*")
end
