-- 3D Transport test with Vacuum and Incident-isotropic BC.
-- SDM: PWLD
-- Test: Max-value=5.28310e-01 and 8.04576e-04
num_procs = 4
if reflecting == nil then
  reflecting = true
end

-- Check num_procs
if check_num_procs == nil and number_of_processes ~= num_procs then
  log.Log(
    LOG_0ERROR,
    "Incorrect amount of processors. "
      .. "Expected "
      .. tostring(num_procs)
      .. ". Pass check_num_procs=false to override if possible."
  )
  os.exit(false)
end

-- Setup mesh
nodes = {}
N = 10
L = 5.0
xmin = -L / 2
dx = L / N
for i = 1, (N + 1) do
  k = i - 1
  nodes[i] = xmin + k * dx
end
znodes = {}
for i = 1, (N / 2 + 1) do
  k = i - 1
  znodes[i] = xmin + k * dx
end

if reflecting then
  meshgen1 = mesh.OrthogonalMeshGenerator.Create({ node_sets = { nodes, nodes, znodes } })
else
  meshgen1 = mesh.OrthogonalMeshGenerator.Create({ node_sets = { nodes, nodes, nodes } })
end
grid = meshgen1:Execute()

-- Set Material IDs
vol0 = logvol.RPPLogicalVolume.Create({ infx = true, infy = true, infz = true })
grid:SetMaterialIDFromLogicalVolume(vol0, 0, true)

num_groups = 21
xs_graphite = xs.LoadFromOpenSn("xs_graphite_pure.xs")

strength = {}
for g = 1, num_groups do
  strength[g] = 0.0
end
mg_src = lbs.VolumetricSource.Create({ block_ids = { 1 }, group_strength = strength })

-- Setup Physics
pquad0 = aquad.CreateProductQuadrature(GAUSS_LEGENDRE_CHEBYSHEV, 2, 2)

lbs_block = {
  mesh = grid,
  num_groups = num_groups,
  groupsets = {
    {
      groups_from_to = { 0, 20 },
      angular_quadrature = pquad0,
      angle_aggregation_type = "single",
      angle_aggregation_num_subsets = 1,
      groupset_num_subsets = 1,
      inner_linear_method = "petsc_gmres",
      l_abs_tol = 1.0e-6,
      l_max_its = 300,
      gmres_restart_interval = 100,
    },
  },
  xs_map = {
    { block_ids = { 0, 1 }, xs = xs_graphite },
  },
}
bsrc = {}
for g = 1, num_groups do
  bsrc[g] = 0.0
end
bsrc[1] = 1.0 / 4.0 / math.pi
lbs_options = {
  boundary_conditions = {
    { name = "xmin", type = "isotropic", group_strength = bsrc },
  },
  scattering_order = 1,
  volumetric_sources = { mg_src },
}
if reflecting then
  table.insert(lbs_options.boundary_conditions, { name = "zmin", type = "reflecting" })
end

phys1 = lbs.DiscreteOrdinatesSolver.Create(lbs_block)
phys1:SetOptions(lbs_options)

-- Initialize and Execute Solver
ss_solver = lbs.SteadyStateSolver.Create({ lbs_solver = phys1 })

ss_solver:Initialize()
ss_solver:Execute()

-- Get field functions
fflist = lbs.GetScalarFieldFunctionList(phys1)

-- Slice plot
--slices = {}
--for k=1,count do
--    slices[k] = fieldfunc.FFInterpolationCreate(SLICE)
--    fieldfunc.SetProperty(slices[k],SLICE_POINT,{x = 0.0, y = 0.0, z = 0.8001})
--    fieldfunc.SetProperty(slices[k],ADD_FIELDFUNCTION,fflist[k])
--    --fieldfunc.SetProperty(slices[k],SLICE_TANGENT,{x = 0.393, y = 1.0-0.393, z = 0})
--    --fieldfunc.SetProperty(slices[k],SLICE_NORMAL,{x = -(1.0-0.393), y = -0.393, z = 0.0})
--    --fieldfunc.SetProperty(slices[k],SLICE_BINORM,{x = 0.0, y = 0.0, z = 1.0})
--    fieldfunc.Initialize(slices[k])
--    fieldfunc.Execute(slices[k])
--    fieldfunc.ExportToPython(slices[k])
--end

-- Volume integrations
ffi1 = fieldfunc.FieldFunctionInterpolationVolume.Create()
curffi = ffi1
curffi:SetOperationType(OP_MAX)
curffi:SetLogicalVolume(vol0)
curffi:AddFieldFunction(fflist[1])

curffi:Initialize()
curffi:Execute()
maxval = curffi:GetValue()

log.Log(LOG_0, string.format("Max-value1=%.5e", maxval))

ffi1 = fieldfunc.FieldFunctionInterpolationVolume.Create()
curffi = ffi1
curffi:SetOperationType(OP_MAX)
curffi:SetLogicalVolume(vol0)
curffi:AddFieldFunction(fflist[20])

curffi:Initialize()
curffi:Execute()
maxval = curffi:GetValue()

log.Log(LOG_0, string.format("Max-value2=%.5e", maxval))

-- Exports
if master_export == nil then
  if reflecting then
    fieldfunc.ExportToVTKMulti(fflist, "ZPhi3DReflected")
  else
    fieldfunc.ExportToVTKMulti(fflist, "ZPhi3D")
  end
end

-- Plots
if location_id == 0 and master_export == nil then
  --os.execute("python ZPFFI00.py")
  ----os.execute("python ZPFFI11.py")
  --local handle = io.popen("python ZPFFI00.py")
  print("Execution completed")
end
