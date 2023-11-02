--############################################### Setup mesh
nodes={}
N=100
L=2.0
xmin = -L/2
dx = L/N
for i=1,(N+1) do
  k=i-1
  nodes[i] = xmin + k*dx
end

meshgen1 = chi_mesh.OrthogonalMeshGenerator.Create({ node_sets = {nodes, nodes} })
chi_mesh.MeshGenerator.Execute(meshgen1)