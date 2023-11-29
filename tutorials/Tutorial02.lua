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

chiRegionExportMeshToVTK(region1,"Mesh")
--############################################### Add material properties


-- Set material properties
PhysicsMaterialAddProperty(material,SCALAR_VALUE,"k")
PhysicsMaterialSetProperty(material,"k",SINGLE_VALUE,1.0)

PhysicsMaterialAddProperty(material,SCALAR_VALUE,"q")
PhysicsMaterialSetProperty(material,"q",SINGLE_VALUE,1.0)


--############################################### Setup Physics
phys1 = DiffusionCreateSolver();
SolverAddRegion(phys1,region1)
SolverSetBasicOption(phys1,"discretization_method","PWLC");
SolverSetBasicOption(phys1,"residual_tolerance",1.0e-6)

--############################################### Initialize and
--                                                Execute Solver
DiffusionInitialize(phys1)
DiffusionExecute(phys1)

----############################################### Visualize the field function
fflist,count = chiGetFieldFunctionList(phys1)
ExportFieldFunctionToVTK(fflist[1],"Tutorial1Output","Temperature")

slice1 = FFInterpolationCreate(SLICE)
FFInterpolationSetProperty(slice1,SLICE_POINT,0.0,0.0,0.0)
FFInterpolationSetProperty(slice1,ADD_FIELDFUNCTION,fflist[1])

FFInterpolationInitialize(slice1)
FFInterpolationExecute(slice1)
FFInterpolationExportPython(slice1)

local handle = io.popen("python ZPFFI00.py")
