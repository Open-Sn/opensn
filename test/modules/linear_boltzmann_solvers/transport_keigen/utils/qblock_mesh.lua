--############################################### Setup mesh
nodes={}
N=40
L=14.0
xmin = 0.0
dx = L/N
for i=1,(N+1) do
    k=i-1
    nodes[i] = xmin + k*dx
end
meshgen1 = chi_mesh.OrthogonalMeshGenerator.Create({ node_sets = {nodes,nodes} })
chi_mesh.MeshGenerator.Execute(meshgen1)

chiVolumeMesherSetMatIDToAll(0)

vol1 = chi_mesh.RPPLogicalVolume.Create
({ xmin=-1000.0,xmax=10.0,ymin=-1000.0,ymax=10.0, infz=true })
chiVolumeMesherSetProperty(MATID_FROMLOGICAL,vol1,1)

