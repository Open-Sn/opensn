-- 1D Diffusion test with Vacuum BCs.
-- SDM: PWLC
-- Test: Max-value=2.50000
num_procs = 1





--############################################### Check num_procs
if (check_num_procs==nil and chi_number_of_processes ~= num_procs) then
    chiLog(LOG_0ERROR,"Incorrect amount of processors. " ..
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
chiVolumeMesherSetMatIDToAll(0)
chiVolumeMesherSetupOrthogonalBoundaries()

--############################################### Add materials
materials = {}
materials[0] = chiPhysicsAddMaterial("Test Material");

chiPhysicsMaterialAddProperty(materials[0],SCALAR_VALUE)
chiPhysicsMaterialSetProperty(materials[0],SCALAR_VALUE,SINGLE_VALUE,1.0)

--############################################### Setup Physics
phys1 = chiDiffusionCreateSolver();
chiSolverSetBasicOption(phys1,"discretization_method","PWLC")
chiSolverSetBasicOption(phys1,"residual_tolerance",1.0e-4)

--############################################### Set boundary conditions
chiDiffusionSetProperty(phys1,"boundary_type","ZMIN","vacuum")
chiDiffusionSetProperty(phys1,"boundary_type","ZMAX","vacuum")

--############################################### Initialize and Execute Solver
chiDiffusionInitialize(phys1)
chiDiffusionExecute(phys1)

--############################################### Get field functions
fftemp,count = chiSolverGetFieldFunctionList(phys1)

--############################################### Line plot
ffi0 = chiFFInterpolationCreate(LINE)
curffi = ffi0;
chiFFInterpolationSetProperty(curffi,LINE_FIRSTPOINT,0.0,0.0,0.0+xmin)
chiFFInterpolationSetProperty(curffi,LINE_SECONDPOINT,0.0,0.0, 2.0+xmin)
chiFFInterpolationSetProperty(curffi,LINE_NUMBEROFPOINTS, 1000)
chiFFInterpolationSetProperty(curffi,ADD_FIELDFUNCTION,fftemp[1])

chiFFInterpolationInitialize(curffi)
chiFFInterpolationExecute(curffi)

--############################################### Volume integrations
vol0 = mesh.RPPLogicalVolume.Create({infx=true, infy=true, infz=true})
ffi1 = chiFFInterpolationCreate(VOLUME)
curffi = ffi1
chiFFInterpolationSetProperty(curffi,OPERATION,OP_MAX)
chiFFInterpolationSetProperty(curffi,LOGICAL_VOLUME,vol0)
chiFFInterpolationSetProperty(curffi,ADD_FIELDFUNCTION,fftemp[1])

chiFFInterpolationInitialize(curffi)
chiFFInterpolationExecute(curffi)
maxval = chiFFInterpolationGetValue(curffi)

chiLog(LOG_0,string.format("Max-value=%.5f", maxval))

--############################################### Exports
if (master_export == nil) then
    chiFFInterpolationExportPython(ffi0)
end

--############################################### Plots
if (location_id == 0 and master_export == nil) then
    local handle = io.popen("python ZLFFI00.py")
end
