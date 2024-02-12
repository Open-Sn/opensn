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
mesh.SetMatIDToAll(0)


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

mesh.SetProperty(BNDRYID_FROMLOGICAL,e_vol,e_bndry)
mesh.SetProperty(BNDRYID_FROMLOGICAL,w_vol,w_bndry)
mesh.SetProperty(BNDRYID_FROMLOGICAL,n_vol,n_bndry)
mesh.SetProperty(BNDRYID_FROMLOGICAL,s_vol,s_bndry)

--############################################### Add material properties
--#### DFEM solver
phys1 = DFEMDiffusionSolverCreate()

DFEMDiffusionSetBCProperty(phys1,"boundary_type",e_bndry,"dirichlet",0.0)
DFEMDiffusionSetBCProperty(phys1,"boundary_type",w_bndry,"dirichlet",0.0)
DFEMDiffusionSetBCProperty(phys1,"boundary_type",n_bndry,"dirichlet",0.0)
DFEMDiffusionSetBCProperty(phys1,"boundary_type",s_bndry,"dirichlet",0.0)

SolverInitialize(phys1)
SolverExecute(phys1)

--############################################### Get field functions
fflist,count = SolverGetFieldFunctionList(phys1)

--############################################### Export VTU
if (master_export == nil) then
    ExportFieldFunctionToVTK(fflist[1],"CFEMDiff2D_analytic_coef","flux")
end

--############################################### Volume integrations
vol0 = mesh.RPPLogicalVolume.Create({infx=true, infy=true, infz=true})

ffvol = FFInterpolationCreate(VOLUME)
FFInterpolationSetProperty(ffvol,OPERATION,OP_MAX)
FFInterpolationSetProperty(ffvol,LOGICAL_VOLUME,vol0)
FFInterpolationSetProperty(ffvol,ADD_FIELDFUNCTION,fflist[1])

FFInterpolationInitialize(ffvol)
FFInterpolationExecute(ffvol)
maxval = FFInterpolationGetValue(ffvol)

Log(LOG_0,string.format("Max-value=%.6f", maxval))
