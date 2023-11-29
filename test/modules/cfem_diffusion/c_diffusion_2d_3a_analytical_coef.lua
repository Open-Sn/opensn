--############################################### Setup mesh
nodes={}
N=100
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
chiVolumeMesherSetMatIDToAll(0)


function D_coef(i,x,y,z)
    return 3.0 + x + y
end
function Q_ext(i,x,y,z)
    return x*x
end
function Sigma_a(i,x,y,z)
    return x*y*y
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

--MeshHandlerExportMeshToVTK("Mesh")

--############################################### Add material properties
--#### CFEM solver
phys1 = chiCFEMDiffusionSolverCreate()

chiCFEMDiffusionSetBCProperty(phys1,"boundary_type",e_bndry,"dirichlet",0.0)
chiCFEMDiffusionSetBCProperty(phys1,"boundary_type",w_bndry,"dirichlet",0.0)
chiCFEMDiffusionSetBCProperty(phys1,"boundary_type",n_bndry,"dirichlet",0.0)
chiCFEMDiffusionSetBCProperty(phys1,"boundary_type",s_bndry,"dirichlet",0.0)

chiSolverInitialize(phys1)
chiSolverExecute(phys1)

--############################################### Get field functions
fflist,count = chiSolverGetFieldFunctionList(phys1)

--############################################### Export VTU
if (master_export == nil) then
    chiExportFieldFunctionToVTK(fflist[1],"CFEMDiff2D_analytic_coef","flux")
end

--############################################### Volume integrations

--############################################### PostProcessors
chi.AggregateNodalValuePostProcessor.Create
({
    name = "maxval",
    field_function = math.floor(fflist[1]),
    operation = "max"
})
chi.ExecutePostProcessors({"maxval"})
