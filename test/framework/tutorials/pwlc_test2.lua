--############################################### Setup mesh
if (nmesh==nil) then nmesh = 10 end

nodes={}
N=nmesh
L=2.0
xmin = -L/2
dx = L/N
for i=1,(N+1) do
    k=i-1
    nodes[i] = xmin + k*dx
end

meshgen1 = mesh.OrthogonalMeshGenerator.Create({ node_sets = {nodes,nodes} })
mesh.MeshGenerator.Execute(meshgen1)

--############################################### Set Material IDs
chiVolumeMesherSetMatIDToAll(0)

function MMS_phi(x,y,z)
    return math.cos(math.pi*x) + math.cos(math.pi*y)
end
function MMS_q(x,y,z)
    return math.pi*math.pi * (math.cos(math.pi*x)+math.cos(math.pi*y))
end

chi_unit_tests.chiSimTest04_PWLC()
chiMPIBarrier()
if (location_id == 0) then
    os.execute("rm CodeTut4_PWLC*")
end
