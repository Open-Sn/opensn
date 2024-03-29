--############################################### Setup mesh
nodes={}
N=32
ds=2.0/N
for i=0,N do
    nodes[i+1] = -1.0 + i*ds
end
meshgen1 = mesh.OrthogonalMeshGenerator.Create({ node_sets = {nodes,nodes,nodes} })
mesh.MeshGenerator.Execute(meshgen1)

material = mat.AddMaterial("Test Material");

-- Set Material IDs
vol0 = logvol.Create(RPP,-1000,1000,-1000,1000,-1000,1000)
mesh.SetMaterialIDFromLogicalVolume(vol0,material)

RegionExportMeshToVTK(region1,"Mesh")
--############################################### Add material properties


-- Set material properties
mat.AddProperty(material, SCALAR_VALUE, "k")
mat.SetProperty(material, "k", SINGLE_VALUE, 1.0)

mat.AddProperty(material, SCALAR_VALUE, "q")
mat.SetProperty(material, "q", SINGLE_VALUE, 1.0)


--############################################### Setup Physics
phys1 = DiffusionCreateSolver();
SolverAddRegion(phys1,region1)
solver.SetBasicOption(phys1,"discretization_method","PWLC");
solver.SetBasicOption(phys1,"residual_tolerance",1.0e-6)

--############################################### Initialize and
--                                                Execute Solver
DiffusionInitialize(phys1)
DiffusionExecute(phys1)

----############################################### Visualize the field function
fflist,count = GetFieldFunctionList(phys1)
fieldfunc.ExportToVTK(fflist[1],"Tutorial1Output","Temperature")

slice1 = fieldfunc.FFInterpolationCreate(SLICE)
fieldfunc.SetProperty(slice1,SLICE_POINT,0.0,0.0,0.0)
fieldfunc.SetProperty(slice1,ADD_FIELDFUNCTION,fflist[1])

fieldfunc.Initialize(slice1)
fieldfunc.Execute(slice1)
fieldfunc.ExportPython(slice1)

local handle = io.popen("python ZPFFI00.py")
