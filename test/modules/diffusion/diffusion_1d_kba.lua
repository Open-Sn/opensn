-- 1D Diffusion test with Vacuum BCs.
-- SDM: PWLC
-- Test: Max-value=2.50000
num_procs = 2
-- KBA-partitioning




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

meshgen1 = mesh.OrthogonalMeshGenerator.Create
({
  node_sets = {nodes},
  partitioner = chi.KBAGraphPartitioner.Create
  ({
    nx = 1, ny=1, nz=2,
    zcuts = {L/2}
  })
})
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
phys1 = DiffusionCreateSolver();
SolverSetBasicOption(phys1,"discretization_method","PWLC")
SolverSetBasicOption(phys1,"residual_tolerance",1.0e-4)

DiffusionSetProperty(phys1,"boundary_type","ZMIN","vacuum")
DiffusionSetProperty(phys1,"boundary_type","ZMAX","vacuum")


--############################################### Initialize and Execute Solver
DiffusionInitialize(phys1)
DiffusionExecute(phys1)

--############################################### Get field functions
fftemp,count = SolverGetFieldFunctionList(phys1)

--############################################### Line plot
ffi0 = FFInterpolationCreate(LINE)
curffi = ffi0;
FFInterpolationSetProperty(curffi,LINE_FIRSTPOINT,0.0,0.0,0.0+xmin)
FFInterpolationSetProperty(curffi,LINE_SECONDPOINT,0.0,0.0, 2.0+xmin)
FFInterpolationSetProperty(curffi,LINE_NUMBEROFPOINTS, 1000)
FFInterpolationSetProperty(curffi,ADD_FIELDFUNCTION,fftemp[1])

FFInterpolationInitialize(curffi)
FFInterpolationExecute(curffi)

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

Log(LOG_0,string.format("Max-value=%.5f", maxval))

--############################################### Exports
if (master_export == nil) then
    FFInterpolationExportPython(ffi0)
end

if (location_id == 0 and master_export == nil) then

    local handle = io.popen("python ZLFFI00.py")
end
