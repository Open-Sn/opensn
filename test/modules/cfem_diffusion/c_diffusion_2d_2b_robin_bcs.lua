--############################################### Setup mesh
nodes={}
N=40
L=2
xmin = -L/2
dx = L/N
for i=1,(N+1) do
    k=i-1
    nodes[i] = xmin + k*dx
end
 
meshgen1 = mesh.OrthogonalMeshGenerator.Create({ node_sets = {nodes,nodes} })
mesh.MeshGenerator.Execute(meshgen1)
 
--############################################### Set Material IDs
vol0 = mesh.RPPLogicalVolume.Create({infx=true, infy=true, infz=true})
VolumeMesherSetProperty(MATID_FROMLOGICAL,vol0,0)

vol1 = mesh.RPPLogicalVolume.Create
({ xmin=-0.5,xmax=0.5,ymin=-0.5,ymax=0.5, infz=true })
VolumeMesherSetProperty(MATID_FROMLOGICAL,vol1,1)

D = {1.0,5.0}
Q = {0.0,1.0}
XSa = {1.0,1.0}
function D_coef(i,x,y,z)
    return D[i+1] -- + x
end
function Q_ext(i,x,y,z)
    return Q[i+1] -- x*x
end
function Sigma_a(i,x,y,z)
    return XSa[i+1]
end

-- Setboundary IDs
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

--############################################### Add material properties
--#### CFEM solver
phys1 = chiCFEMDiffusionSolverCreate()

chiSolverSetBasicOption(phys1, "residual_tolerance", 1E-8)

chiCFEMDiffusionSetBCProperty(phys1,"boundary_type",e_bndry,"robin", 0.25, 0.5, 0.0)
chiCFEMDiffusionSetBCProperty(phys1,"boundary_type",n_bndry,"reflecting")
chiCFEMDiffusionSetBCProperty(phys1,"boundary_type",s_bndry,"reflecting")
chiCFEMDiffusionSetBCProperty(phys1,"boundary_type",w_bndry,"robin", 0.25, 0.5, 0.1)

SolverInitialize(phys1)
chiSolverExecute(phys1)

--############################################### Get field functions
fflist,count = chiSolverGetFieldFunctionList(phys1)

--############################################### Export VTU
if (master_export == nil) then
    chiExportFieldFunctionToVTK(fflist[1],"CFEMDiff2D_RobinRefl","flux")
end

--############################################### Volume integrations

--############################################### PostProcessors
chi.CellVolumeIntegralPostProcessor.Create
({
    name = "avgval",
    field_function = math.floor(fflist[1]),
    compute_volume_average = true
})
chi.ExecutePostProcessors({"avgval"})
