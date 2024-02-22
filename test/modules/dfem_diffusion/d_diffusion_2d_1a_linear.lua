--############################################### Setup mesh
nodes={}
N=10
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
mesh.SetUniformMaterialID(0)

D = {1.0}
Q = {0.0}
XSa = {0.0}
function D_coef(i,x,y,z)
    return D[i+1]
end
function Q_ext(i,x,y,z)
    return Q[i+1]
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

mesh.SetBoundaryIDFromLogicalVolume(e_vol,e_bndry)
mesh.SetBoundaryIDFromLogicalVolume(w_vol,w_bndry)
mesh.SetBoundaryIDFromLogicalVolume(n_vol,n_bndry)
mesh.SetBoundaryIDFromLogicalVolume(s_vol,s_bndry)

--############################################### Add material properties
--#### DFEM solver
phys1 = DFEMDiffusionSolverCreate()

solver.SetBasicOption(phys1, "residual_tolerance", 1E-8)

DFEMDiffusionSetBCProperty(phys1,"boundary_type",e_bndry,"robin", 0.25, 0.5, 0.0)
DFEMDiffusionSetBCProperty(phys1,"boundary_type",n_bndry,"reflecting")
DFEMDiffusionSetBCProperty(phys1,"boundary_type",s_bndry,"reflecting")
DFEMDiffusionSetBCProperty(phys1,"boundary_type",w_bndry,"robin", 0.25, 0.5, 1.0)

solver.Initialize(phys1)
solver.Execute(phys1)


--############################################### Get field functions
fflist,count = solver.GetFieldFunctionList(phys1)

--############################################### Export VTU
if (master_export == nil) then
    ExportFieldFunctionToVTK(fflist[1],"DFEMDiff2D_linear")
end

--############################################### Line plot
cline = FFInterpolationCreate(LINE)
FFInterpolationSetProperty(cline,LINE_FIRSTPOINT,-L/2, 0.0, 0.0)
FFInterpolationSetProperty(cline,LINE_SECONDPOINT,L/2, 0.0, 0.0)
FFInterpolationSetProperty(cline,LINE_NUMBEROFPOINTS, 50)
FFInterpolationSetProperty(cline,ADD_FIELDFUNCTION,fflist[1])

FFInterpolationInitialize(cline)
FFInterpolationExecute(cline)

if (master_export == nil) then
    FFInterpolationExportPython(cline)
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

log.Log(LOG_0,string.format("Max-value=%.6f", maxval))
