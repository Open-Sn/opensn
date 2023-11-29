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
vol0 = chiLogicalVolumeCreate(RPP,-1000,1000,-1000,1000,-1000,1000)
VolumeMesherSetProperty(MATID_FROMLOGICAL,vol0,material)
-- Setboundary IDs
-- xmin,xmax,ymin,ymax,zmin,zmax
e_vol = chiLogicalVolumeCreate(RPP,0.99999,1000,-1000,1000,-1000,1000)
w_vol = chiLogicalVolumeCreate(RPP,-1000,0.00001,-1000,1000,-1000,1000)
n_vol = chiLogicalVolumeCreate(RPP,-1000,1000,0.99999,1000,-1000,1000)
s_vol = chiLogicalVolumeCreate(RPP,-1000,1000,-1000,0.00001,-1000,1000)

e_bndry = 0
w_bndry = 1
n_bndry = 2
s_bndry = 3

VolumeMesherSetProperty(BNDRYID_FROMLOGICAL,e_vol,e_bndry)
VolumeMesherSetProperty(BNDRYID_FROMLOGICAL,w_vol,w_bndry)
VolumeMesherSetProperty(BNDRYID_FROMLOGICAL,n_vol,n_bndry)
VolumeMesherSetProperty(BNDRYID_FROMLOGICAL,s_vol,s_bndry)

MeshHandlerExportMeshToVTK("Mesh")

--############################################### Add material properties
-- Set material properties
PhysicsMaterialAddProperty(material,SCALAR_VALUE,"D")
chiPhysicsMaterialSetProperty(material,"D",SINGLE_VALUE,1.0)

PhysicsMaterialAddProperty(material,SCALAR_VALUE,"q")
chiPhysicsMaterialSetProperty(material,"q",SINGLE_VALUE,0.0)

--############################################### Setup Physics
phys1 = chiDiffusionCreateSolver()
chiSolverSetBasicOption(phys1,"discretization_method","PWLD_MIP");
chiSolverSetBasicOption(phys1,"residual_tolerance",1.0e-6)
chiDiffusionSetProperty(phys1,"boundary_type",e_bndry,"robin", 0.25, 0.5, 0.0)
chiDiffusionSetProperty(phys1,"boundary_type",n_bndry,"reflecting")
chiDiffusionSetProperty(phys1,"boundary_type",s_bndry,"reflecting")
chiDiffusionSetProperty(phys1,"boundary_type",w_bndry,"robin", 0.25, 0.5, 1.0)

--############################################### Initialize and
--                                                Execute Solver
chiDiffusionInitialize(phys1)
chiDiffusionExecute(phys1)

----############################################### Visualize the field function
fflist,count = chiGetFieldFunctionList(phys1)
chiExportFieldFunctionToVTK(fflist[1],"Tutorial4_Diff_PWLD_Output","Flux_Diff_IP")
