--############################################### Setup mesh
if (nmesh==nil) then nmesh = 10 end

nodes={}
N=nmesh
L=2.0
xmin = -L/2
--xmin = 0.0
dx = L/N
for i=1,(N+1) do
    k=i-1
    nodes[i] = xmin + k*dx
end

meshgen1 = mesh.OrthogonalMeshGenerator.Create({ node_sets = {nodes,nodes} })
mesh.MeshGenerator.Execute(meshgen1)

--############################################### Set Material IDs
mesh.SetUniformMaterialID(0)

mesh.SetupOrthogonalBoundaries()

function MMS_phi(x,y,z)
    return math.cos(math.pi*x) + math.cos(math.pi*y)
end
function MMS_q(x,y,z)
    return math.pi*math.pi * (math.cos(math.pi*x)+math.cos(math.pi*y))
end

unit_tests.acceleration_Diffusion_CFEM();
MPIBarrier()
if (location_id == 0) then
    --os.execute("rm SimTest_92*")
end

--[0]  Iteration     0   1.000e+00
--[0]  Iteration     1   2.016e+02
--[0]  Iteration     2   1.941e+00
--[0]  Iteration     3   1.294e+00
--[0]  Iteration     4   3.890e-01
--[0]  Iteration     5   2.887e-02
--[0]  Iteration     6   1.239e-03
--[0]  Iteration     7   4.076e-05
--[0]  Iteration     8   1.119e-06
--[0]  Iteration     9   2.955e-08
