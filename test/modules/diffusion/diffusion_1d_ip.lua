-- 1D Diffusion test with Dirichlet BCs
-- SMD: PWLD
-- Test: Max-value=0.5006523132
num_procs = 2





--############################################### Check num_procs
if (check_num_procs==nil and number_of_processes ~= num_procs) then
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
phys1 = chiDiffusionCreateSolver()
chiSolverSetBasicOption(phys1,"discretization_method","PWLD_MIP")
chiSolverSetBasicOption(phys1,"residual_tolerance",1.0e-6)

--############################################### Set boundary conditions
chiDiffusionSetProperty(phys1,"boundary_type","ZMIN","dirichlet",0.0)
chiDiffusionSetProperty(phys1,"boundary_type","ZMAX","dirichlet",0.0)

--############################################### Initialize and Execute Solver
chiDiffusionInitialize(phys1)
chiDiffusionExecute(phys1)

--############################################### Get field functions
fftemp,count = chiSolverGetFieldFunctionList(phys1)

--############################################### Line plot
line0 = chiFFInterpolationCreate(LINE)
chiFFInterpolationSetProperty(line0,LINE_FIRSTPOINT,0.1,0.0,0.0)
chiFFInterpolationSetProperty(line0,LINE_SECONDPOINT,0.1,0.0, 2.0)
chiFFInterpolationSetProperty(line0,LINE_NUMBEROFPOINTS, 100)
chiFFInterpolationSetProperty(line0,ADD_FIELDFUNCTION,fftemp[1])

chiFFInterpolationInitialize(line0)
chiFFInterpolationExecute(line0)

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

chiLog(LOG_0,string.format("Max-value=%.10f", maxval))

--############################################### Exports
if (master_export == nil) then
    chiFFInterpolationExportPython(line0)
end

--############################################### Plots
if ((location_id == 0) and (master_export == nil)) then
    local handle = io.popen("python ZLFFI00.py")
end
