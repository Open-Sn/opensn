--############################################### Setup mesh
nodes={}
N=40
L=1
xmin = 0
dx = L/N
for i=1,(N+1) do
    k=i-1
    nodes[i] = xmin + k*dx
end
 
meshgen1 = mesh.OrthogonalMeshGenerator.Create({ node_sets = {nodes,nodes} })
mesh.MeshGenerator.Execute(meshgen1)
 
--############################################### Set Material IDs
VolumeMesherSetMatIDToAll(0)

-- governing law: -(u_xx + u_yy) = q, on domain [0,1]x[0,1]
-- when the exact solution is chosen u(x,y) = sin(pi.x) * sin(pi.y)
-- this automatically gives:
--    boundary = zero-Dirichlet on all 4 sides
--    volumetric source term: q(,x) = 2*pi*pi * sin(pi.x) * sin(pi.y)
-- the factor 2 is the dim of the problem
function D_coef(i,x,y,z)
    return 1.0
end
function Q_ext(i,x,y,z)
    return 2.*math.pi*math.pi * math.sin(math.pi*x) * math.sin(math.pi*y)
end
function Sigma_a(i,x,y,z)
    return 0.0
end
function MMS_phi(i,x,y,z)
    return math.sin(math.pi*x) * math.sin(math.pi*y)
end

-- Set boundary IDs
-- xmin,xmax,ymin,ymax,zmin,zmax
e_vol = mesh.RPPLogicalVolume.Create({xmin=0.99999,xmax=1000.0  , infy=true, infz=true})
w_vol = mesh.RPPLogicalVolume.Create({xmin=-1000.0,xmax=-0.99999, infy=true, infz=true})
n_vol = mesh.RPPLogicalVolume.Create({ymin=0.99999,ymax=1000.0  , infx=true, infz=true})
s_vol = mesh.RPPLogicalVolume.Create({ymin=-1000.0,ymax=-0.99999, infx=true, infz=true})

e_bndry = 0
w_bndry = 1
n_bndry = 2
s_bndry = 3

VolumeMesherSetProperty(BNDRYID_FROMLOGICAL,e_vol,e_bndry)
VolumeMesherSetProperty(BNDRYID_FROMLOGICAL,w_vol,w_bndry)
VolumeMesherSetProperty(BNDRYID_FROMLOGICAL,n_vol,n_bndry)
VolumeMesherSetProperty(BNDRYID_FROMLOGICAL,s_vol,s_bndry)

--############################################### Call Lua Sim Test
SimTest_IP_MMS_L2error() --simtest_IP_MMS_L2_handle becomes available here

--############################################### Export VTU
if (master_export == nil) then
    chiExportFieldFunctionToVTK(simtest_IP_MMS_L2_handle,"DFEMDiff2D_MMS","flux")
end

----############################################### Volume integrations
vol0 = mesh.RPPLogicalVolume.Create({infx=true, infy=true, infz=true})
--
ffvol = chiFFInterpolationCreate(VOLUME)
chiFFInterpolationSetProperty(ffvol,OPERATION,OP_MAX)
chiFFInterpolationSetProperty(ffvol,LOGICAL_VOLUME,vol0)
chiFFInterpolationSetProperty(ffvol,ADD_FIELDFUNCTION,simtest_IP_MMS_L2_handle)
--
chiFFInterpolationInitialize(ffvol)
chiFFInterpolationExecute(ffvol)
maxval = chiFFInterpolationGetValue(ffvol)

chiLog(LOG_0,string.format("Max-value=%.6f", maxval))
