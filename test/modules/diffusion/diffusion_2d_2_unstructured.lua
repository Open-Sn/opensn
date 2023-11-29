-- 2D Diffusion test with Dirichlet BCs.
-- SDM: PWLC
-- Test: Max-value=0.30384
num_procs = 4





--############################################### Check num_procs
if (check_num_procs==nil and number_of_processes ~= num_procs) then
    Log(LOG_0ERROR,"Incorrect amount of processors. " ..
                      "Expected "..tostring(num_procs)..
                      ". Pass check_num_procs=false to override if possible.")
    os.exit(false)
end

--############################################### Setup mesh
meshgen1 = mesh.MeshGenerator.Create
({
  inputs =
  {
    mesh.FromFileMeshGenerator.Create
    ({
      filename = "../../../resources/TestMeshes/TriangleMesh2x2.obj"
    })
  }
})
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
phys1 = DiffusionCreateSolver()
SolverSetBasicOption(phys1,"discretization_method","PWLC")
SolverSetBasicOption(phys1,"residual_tolerance",1.0e-6)

--############################################### Set boundary conditions
--DiffusionSetProperty(phys1,"boundary_type",0,"reflecting",1.0)
--DiffusionSetProperty(phys1,"boundary_type",1,"vacuum",2.0)
--DiffusionSetProperty(phys1,"boundary_type",2,"reflecting",3.0)
--DiffusionSetProperty(phys1,"boundary_type",3,"vacuum",4.0)

--############################################### Initialize and Execute Solver
DiffusionInitialize(phys1)
DiffusionExecute(phys1)

--############################################### Get field functions
fftemp,count = SolverGetFieldFunctionList(phys1)

--############################################### Slice plot
slice2 = FFInterpolationCreate(SLICE)
FFInterpolationSetProperty(slice2,SLICE_POINT,0.0,0.0,0.025)
FFInterpolationSetProperty(slice2,ADD_FIELDFUNCTION,fftemp[1])

FFInterpolationInitialize(slice2)
FFInterpolationExecute(slice2)

--############################################### Line plot
line0 = FFInterpolationCreate(LINE)
FFInterpolationSetProperty(line0,LINE_FIRSTPOINT,-1.0,0.01,0.0)
FFInterpolationSetProperty(line0,LINE_SECONDPOINT, 1.0,0.01,0.0)
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

Log(LOG_0,string.format("Max-value=%.5f", maxval))

--############################################### Exports
if (master_export == nil) then
    FFInterpolationExportPython(slice2)
    FFInterpolationExportPython(line0)

    ExportFieldFunctionToVTK(fftemp,"ZPhi")
end

--############################################### Plots
if ((master_export == nil) and (location_id == 0)) then
    local handle = io.popen("python3 ZPFFI00.py")
    local handle = io.popen("python3 ZLFFI10.py")
end
