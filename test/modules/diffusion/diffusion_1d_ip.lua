-- 1D Diffusion test with Dirichlet BCs
-- SMD: PWLD
-- Test: Max-value=0.5006523132
num_procs = 2





--############################################### Check num_procs
if (check_num_procs==nil and number_of_processes ~= num_procs) then
    Log(LOG_0ERROR,"Incorrect amount of processors. " ..
                      "Expected "..tostring(num_procs)..
                      ". Pass check_num_procs=false to override if possible.")
    os.exit(false)
end

--############################################### Setup mesh
nodes={}
N=100
L=2.0
xmin = 0.0
dx = L/N
for i=1,(N+1) do
    k=i-1
    nodes[i] = xmin + k*dx
end

meshgen1 = mesh.OrthogonalMeshGenerator.Create({ node_sets = {nodes} })
mesh.MeshGenerator.Execute(meshgen1)

--############################################### Set Material IDs
VolumeMesherSetMatIDToAll(0)
VolumeMesherSetupOrthogonalBoundaries()

--############################################### Add materials
materials = {}
materials[0] = PhysicsAddMaterial("Test Material");

PhysicsMaterialAddProperty(materials[0],SCALAR_VALUE)
PhysicsMaterialSetProperty(materials[0],SCALAR_VALUE,SINGLE_VALUE,1.0)

--############################################### Setup Physics
phys1 = DiffusionCreateSolver()
SolverSetBasicOption(phys1,"discretization_method","PWLD_MIP")
SolverSetBasicOption(phys1,"residual_tolerance",1.0e-6)

--############################################### Set boundary conditions
DiffusionSetProperty(phys1,"boundary_type","ZMIN","dirichlet",0.0)
DiffusionSetProperty(phys1,"boundary_type","ZMAX","dirichlet",0.0)

--############################################### Initialize and Execute Solver
DiffusionInitialize(phys1)
DiffusionExecute(phys1)

--############################################### Get field functions
fftemp,count = SolverGetFieldFunctionList(phys1)

--############################################### Line plot
line0 = FFInterpolationCreate(LINE)
FFInterpolationSetProperty(line0,LINE_FIRSTPOINT,0.1,0.0,0.0)
FFInterpolationSetProperty(line0,LINE_SECONDPOINT,0.1,0.0, 2.0)
FFInterpolationSetProperty(line0,LINE_NUMBEROFPOINTS, 100)
FFInterpolationSetProperty(line0,ADD_FIELDFUNCTION,fftemp[1])

FFInterpolationInitialize(line0)
FFInterpolationExecute(line0)

--############################################### Volume integrations
vol0 = mesh.RPPLogicalVolume.Create({infx=true, infy=true, infz=true})
ffi1 = FFInterpolationCreate(VOLUME)
curffi = ffi1
FFInterpolationSetProperty(curffi,OPERATION,OP_MAX)
FFInterpolationSetProperty(curffi,LOGICAL_VOLUME,vol0)
FFInterpolationSetProperty(curffi,ADD_FIELDFUNCTION,fftemp[1])

FFInterpolationInitialize(curffi)
FFInterpolationExecute(curffi)
maxval = FFInterpolationGetValue(curffi)

Log(LOG_0,string.format("Max-value=%.10f", maxval))

--############################################### Exports
if (master_export == nil) then
    FFInterpolationExportPython(line0)
end

--############################################### Plots
if ((location_id == 0) and (master_export == nil)) then
    local handle = io.popen("python ZLFFI00.py")
end
