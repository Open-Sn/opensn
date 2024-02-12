-- 3D Transport test with Vacuum and Incident-isotropic BC.
-- SDM: PWLD
-- Test: Max-value=5.28310e-01 and 8.04576e-04
num_procs = 4
if (reflecting == nil) then reflecting = true end




--############################################### Check num_procs
if (check_num_procs==nil and number_of_processes ~= num_procs) then
  Log(LOG_0ERROR,"Incorrect amount of processors. " ..
    "Expected "..tostring(num_procs)..
    ". Pass check_num_procs=false to override if possible.")
  os.exit(false)
end

--############################################### Setup mesh
nodes={}
N=10
L=5
xmin = -L/2
dx = L/N
for i=1,(N+1) do
  k=i-1
  nodes[i] = xmin + k*dx
end
znodes={}
for i=1,(N/2+1) do
  k=i-1
  znodes[i] = xmin + k*dx
end
if (reflecting) then
  meshgen1 = mesh.OrthogonalMeshGenerator.Create({ node_sets = {nodes,nodes,znodes} })
else
  meshgen1 = mesh.OrthogonalMeshGenerator.Create({ node_sets = {nodes,nodes,nodes} })
end
mesh.MeshGenerator.Execute(meshgen1)

--############################################### Set Material IDs
vol0 = mesh.RPPLogicalVolume.Create({infx=true, infy=true, infz=true})
mesh.SetMatIDToAll(0)

--############################################### Add materials
materials = {}
materials[1] = PhysicsAddMaterial("Test Material");

PhysicsMaterialAddProperty(materials[1],TRANSPORT_XSECTIONS)

PhysicsMaterialAddProperty(materials[1],ISOTROPIC_MG_SOURCE)


num_groups = 21
PhysicsMaterialSetProperty(materials[1],TRANSPORT_XSECTIONS,
  OPENSN_XSFILE,"../transport_steady/xs_graphite_pure.xs")

src={}
for g=1,num_groups do
  src[g] = 0.0
end
src[1] = 1.0
PhysicsMaterialSetProperty(materials[1],ISOTROPIC_MG_SOURCE,FROM_ARRAY,src)

--############################################### Setup Physics
pquad0 = CreateProductQuadrature(GAUSS_LEGENDRE_CHEBYSHEV,2, 2)

lbs_block =
{
  num_groups = num_groups,
  groupsets =
  {
    {
      groups_from_to = {0, 20},
      angular_quadrature_handle = pquad0,
      angle_aggregation_type = "single",
      angle_aggregation_num_subsets = 1,
      groupset_num_subsets = 1,
      inner_linear_method = "gmres",
      l_abs_tol = 1.0e-6,
      l_max_its = 300,
      gmres_restart_interval = 100,
    },
  }
}

lbs_options =
{
  scattering_order = 1,
}
if (reflecting) then
  lbs_options.boundary_conditions = {{name = "zmax", type = "reflecting"}}
end

phys1 = lbs.DiffusionDFEMSolver.Create(lbs_block)
lbs.SetOptions(phys1, lbs_options)

--############################################### Initialize and Execute Solver
ss_solver = lbs.SteadyStateSolver.Create({lbs_solver_handle = phys1})

SolverInitialize(ss_solver)
SolverExecute(ss_solver)

--############################################### Get field functions
fflist,count = LBSGetScalarFieldFunctionList(phys1)

--############################################### Slice plot
--slices = {}
--for k=1,count do
--    slices[k] = FFInterpolationCreate(SLICE)
--    FFInterpolationSetProperty(slices[k],SLICE_POINT,0.0,0.0,0.8001)
--    FFInterpolationSetProperty(slices[k],ADD_FIELDFUNCTION,fflist[k])
--    --FFInterpolationSetProperty(slices[k],SLICE_TANGENT,0.393,1.0-0.393,0)
--    --FFInterpolationSetProperty(slices[k],SLICE_NORMAL,-(1.0-0.393),-0.393,0.0)
--    --FFInterpolationSetProperty(slices[k],SLICE_BINORM,0.0,0.0,1.0)
--    FFInterpolationInitialize(slices[k])
--    FFInterpolationExecute(slices[k])
--    FFInterpolationExportPython(slices[k])
--end

--############################################### Volume integrations
ffi1 = FFInterpolationCreate(VOLUME)
curffi = ffi1
FFInterpolationSetProperty(curffi,OPERATION,OP_MAX)
FFInterpolationSetProperty(curffi,LOGICAL_VOLUME,vol0)
FFInterpolationSetProperty(curffi,ADD_FIELDFUNCTION,fflist[1])

FFInterpolationInitialize(curffi)
FFInterpolationExecute(curffi)
maxval = FFInterpolationGetValue(curffi)

Log(LOG_0,string.format("Max-value1=%.5e", maxval))

ffi1 = FFInterpolationCreate(VOLUME)
curffi = ffi1
FFInterpolationSetProperty(curffi,OPERATION,OP_MAX)
FFInterpolationSetProperty(curffi,LOGICAL_VOLUME,vol0)
FFInterpolationSetProperty(curffi,ADD_FIELDFUNCTION,fflist[20])

FFInterpolationInitialize(curffi)
FFInterpolationExecute(curffi)
maxval = FFInterpolationGetValue(curffi)

Log(LOG_0,string.format("Max-value2=%.5e", maxval))

--############################################### Exports
if (master_export == nil) then
  if (reflecting) then
    ExportMultiFieldFunctionToVTK(fflist,"ZPhi3DReflected")
  else
    ExportMultiFieldFunctionToVTK(fflist,"ZPhi3D")
  end
end

--############################################### Plots
if (location_id == 0 and master_export == nil) then

  --os.execute("python ZPFFI00.py")
  ----os.execute("python ZPFFI11.py")
  --local handle = io.popen("python ZPFFI00.py")
  print("Execution completed")
end

-- DO
--[0]  Max-value1=2.52092e+00
--[0]  Max-value2=5.79100e-03

-- MGDiffusion
--[0]  Max-value1=2.74873e-01
--[0]  Max-value2=9.47508e-05
