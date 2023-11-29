-- 3D Diffusion test with Dirichlet and Reflecting BCs.
-- SDM: PWLC
-- Test: Max-value=0.29480
num_procs = 1





--############################################### Check num_procs
if (check_num_procs==nil and number_of_processes ~= num_procs) then
    chiLog(LOG_0ERROR,"Incorrect amount of processors. " ..
                      "Expected "..tostring(num_procs)..
                      ". Pass check_num_procs=false to override if possible.")
    os.exit(false)
end

--############################################### Setup mesh
nodes={}
N=32
L=2.0
xmin = -1.0
dx = L/N
for i=1,(N+1) do
    k=i-1
    nodes[i] = xmin + k*dx
end
znodes={}
for i=1,11 do
    znodes[i] = 0.0 + (i-1)*0.2/10
end

meshgen1 = mesh.OrthogonalMeshGenerator.Create({ node_sets = {nodes,nodes,znodes} })
mesh.MeshGenerator.Execute(meshgen1)

--############################################### Set Material IDs
vol0 = mesh.RPPLogicalVolume.Create({infx=true, infy=true, infz=true})
VolumeMesherSetProperty(MATID_FROMLOGICAL,vol0,0)
VolumeMesherSetupOrthogonalBoundaries()

--############################################### Add materials
materials = {}
materials[0] = PhysicsAddMaterial("Test Material");

PhysicsMaterialAddProperty(materials[0],SCALAR_VALUE)
PhysicsMaterialSetProperty(materials[0],SCALAR_VALUE,SINGLE_VALUE,1.0)

--############################################### Setup Physics
phys1 = chiDiffusionCreateSolver();
SolverSetBasicOption(phys1,"discretization_method","PWLC")
SolverSetBasicOption(phys1,"residual_tolerance",1.0e-6)

--############################################### Set boundary conditions
chiDiffusionSetProperty(phys1,"boundary_type","ZMIN","reflecting")
chiDiffusionSetProperty(phys1,"boundary_type","ZMAX","reflecting")

--############################################### Initialize and Execute Solver
chiDiffusionInitialize(phys1)
chiDiffusionExecute(phys1)

--############################################### Get field functions
fftemp,count = SolverGetFieldFunctionList(phys1)

--############################################### Slice plot
slice1 = FFInterpolationCreate(SLICE)
--FFInterpolationSetProperty(slice1,SLICE_POINT,0.0,0.0,0.025)
FFInterpolationSetProperty(slice1,SLICE_POINT,0.008,0.0,0.0)
FFInterpolationSetProperty(slice1,SLICE_BINORM,0.0,0.0,1.0)
FFInterpolationSetProperty(slice1,SLICE_TANGENT,0.0,-1.0,0.0)
FFInterpolationSetProperty(slice1,SLICE_NORMAL,1.0,0.0,0.0)
FFInterpolationSetProperty(slice1,ADD_FIELDFUNCTION,fftemp[1])

FFInterpolationInitialize(slice1)
FFInterpolationExecute(slice1)

slice2 = FFInterpolationCreate(SLICE)
FFInterpolationSetProperty(slice2,SLICE_POINT,0.0,0.0,0.025)
FFInterpolationSetProperty(slice2,ADD_FIELDFUNCTION,fftemp[1])

FFInterpolationInitialize(slice2)
FFInterpolationExecute(slice2)

--############################################### Line plot
line0 = FFInterpolationCreate(LINE)
FFInterpolationSetProperty(line0,LINE_FIRSTPOINT,-1.0,0.0,0.025)
FFInterpolationSetProperty(line0,LINE_SECONDPOINT, 1.0,0.0,0.025)
FFInterpolationSetProperty(line0,LINE_NUMBEROFPOINTS, 100)
FFInterpolationSetProperty(line0,ADD_FIELDFUNCTION,fftemp[1])

FFInterpolationInitialize(line0)
FFInterpolationExecute(line0)

--############################################### Volume integrations
ffi1 = FFInterpolationCreate(VOLUME)
curffi = ffi1
FFInterpolationSetProperty(curffi,OPERATION,OP_MAX)
FFInterpolationSetProperty(curffi,LOGICAL_VOLUME,vol0)
FFInterpolationSetProperty(curffi,ADD_FIELDFUNCTION,fftemp[1])

FFInterpolationInitialize(curffi)
FFInterpolationExecute(curffi)
maxval = FFInterpolationGetValue(curffi)

chiLog(LOG_0,string.format("Max-value=%.5f", maxval))

--############################################### Exports
if (master_export == nil) then
    FFInterpolationExportPython(slice1)
    FFInterpolationExportPython(slice2)
    FFInterpolationExportPython(line0)
    chiExportFieldFunctionToVTK(fftemp,"ZPhi")
end

--############################################### Plots
if (location_id == 0 and master_export == nil) then
    local handle = io.popen("python3 ZPFFI00.py")
    local handle = io.popen("python3 ZPFFI10.py")
    local handle = io.popen("python3 ZLFFI20.py")
    print("Execution completed")
end
