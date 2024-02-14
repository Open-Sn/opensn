-- 2D Diffusion test with Vacuum and Robin BCs.
-- SDM: PWLD


--############################################### Setup mesh
nodes={}
N=32
L=1.0
ds=L/N
xmin = 0.0
for i=0,N do
    nodes[i+1] = xmin + i*ds
end

meshgen1 = mesh.OrthogonalMeshGenerator.Create({ node_sets = {nodes,nodes} })
mesh.MeshGenerator.Execute(meshgen1)

--############################################### Set Material IDs
material = PhysicsAddMaterial("Homogenous_Material");
-- Set Material IDs
vol0 = LogicalVolumeCreate(RPP,-1000,1000,-1000,1000,-1000,1000)
mesh.SetMaterialIDFromLogicalVolume(vol0,material)
-- Setboundary IDs
-- xmin,xmax,ymin,ymax,zmin,zmax
e_vol = LogicalVolumeCreate(RPP,0.99999,1000,-1000,1000,-1000,1000)
w_vol = LogicalVolumeCreate(RPP,-1000,0.00001,-1000,1000,-1000,1000)
n_vol = LogicalVolumeCreate(RPP,-1000,1000,0.99999,1000,-1000,1000)
s_vol = LogicalVolumeCreate(RPP,-1000,1000,-1000,0.00001,-1000,1000)

e_bndry = 0
w_bndry = 1
n_bndry = 2
s_bndry = 3

mesh.SetBoundaryIDFromLogicalVolume(e_vol,e_bndry)
mesh.SetBoundaryIDFromLogicalVolume(w_vol,w_bndry)
mesh.SetBoundaryIDFromLogicalVolume(n_vol,n_bndry)
mesh.SetBoundaryIDFromLogicalVolume(s_vol,s_bndry)

mesh.ExportToVTK("Mesh")

--############################################### Add material properties
-- Set material properties
PhysicsMaterialAddProperty(material,SCALAR_VALUE,"D")
PhysicsMaterialSetProperty(material,"D",SINGLE_VALUE,1.0)

PhysicsMaterialAddProperty(material,SCALAR_VALUE,"q")
PhysicsMaterialSetProperty(material,"q",SINGLE_VALUE,0.0)

--############################################### Setup Physics
phys1 = DiffusionCreateSolver()
SolverSetBasicOption(phys1,"discretization_method","PWLD_MIP");
SolverSetBasicOption(phys1,"residual_tolerance",1.0e-6)
DiffusionSetProperty(phys1,"boundary_type",e_bndry,"robin", 0.25, 0.5, 0.0)
DiffusionSetProperty(phys1,"boundary_type",n_bndry,"reflecting")
DiffusionSetProperty(phys1,"boundary_type",s_bndry,"reflecting")
DiffusionSetProperty(phys1,"boundary_type",w_bndry,"robin", 0.25, 0.5, 1.0)

--############################################### Initialize and
--                                                Execute Solver
DiffusionInitialize(phys1)
DiffusionExecute(phys1)

----############################################### Visualize the field function
fflist,count = GetFieldFunctionList(phys1)
ExportFieldFunctionToVTK(fflist[1],"Tutorial4_Diff_PWLD_Output","Flux_Diff_IP")
