--############################################### Setup mesh
nodes={}
N=32
ds=2.0/N
for i=0,N do
    nodes[i+1] = -1.0 + i*ds
end
meshgen1 = mesh.OrthogonalMeshGenerator.Create({ node_sets = {nodes,nodes,nodes} })
mesh.MeshGenerator.Execute(meshgen1)

material = PhysicsAddMaterial("Test Material");

-- Set Material IDs
vol0 = chiLogicalVolumeCreate(RPP,-1000,1000,-1000,1000,-1000,1000)
VolumeMesherSetProperty(MATID_FROMLOGICAL,vol0,material)

VolumeMesherSetupOrthogonalBoundaries()

MeshHandlerExportMeshToVTK("Mesh")
--############################################### Add material properties


-- Set material properties
PhysicsMaterialAddProperty(material,SCALAR_VALUE,"k")
PhysicsMaterialSetProperty(material,"k",SINGLE_VALUE,1.0)

PhysicsMaterialAddProperty(material,SCALAR_VALUE,"q")
PhysicsMaterialSetProperty(material,"q",SINGLE_VALUE,1.0)


--############################################### Setup Physics
phys1 = DiffusionCreateSolver()
SolverSetBasicOption(phys1,"discretization_method","PWLC");
SolverSetBasicOption(phys1,"residual_tolerance",1.0e-6)
DiffusionSetProperty(phys1,"boundary_type",4,"reflecting")
DiffusionSetProperty(phys1,"boundary_type",5,"reflecting")

--############################################### Initialize and
--                                                Execute Solver
DiffusionInitialize(phys1)
DiffusionExecute(phys1)

----############################################### Visualize the field function
fflist,count = GetFieldFunctionList(phys1)
ExportFieldFunctionToVTK(fflist[1],"Tutorial1Output","Temperature")
