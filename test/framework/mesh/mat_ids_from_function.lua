-- Setup mesh
nodes={}
N=10
L=5
xmin = -L/2
dx = L/N
for i=1,(N+1) do
    k=i-1
    nodes[i] = xmin + k*dx
end

meshgen1 = mesh.OrthogonalMeshGenerator.Create({ node_sets = {nodes,nodes,nodes} })
mesh.MeshGenerator.Execute(meshgen1)

--Sets a middle square to material 1
function mat_id(x,y,z,cur_id)
    if (math.abs(x)<L/10 and math.abs(y)<L/10) then
        return 1
    end
    return cur_id
end

mesh.SetMaterialIDFromFunction("mat_id")

mesh.ExportToVTK("new_mat_ids")
