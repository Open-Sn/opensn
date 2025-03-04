-- 1D Transport test with Vacuum and Incident-isotropic BC.
-- SDM: PWLD
-- Test: Max-value=0.49903 and 7.18243e-4
num_procs = 3

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
N = 100
L = 30.0
xmin = 0.0
dx = L / N
for i = 1, (N + 1) do
  k = i - 1
  nodes[i] = xmin + k * dx
end

meshgen1 = mesh.OrthogonalMeshGenerator.Create({ node_sets = { nodes } })
grid = meshgen1:Execute()

-- Set Material IDs
grid:SetUniformMaterialID(0)

num_groups = 168
xs_3_170 = xs.LoadFromOpenSn("xs_3_170.xs")

strength = {}
for g = 1, num_groups do
  strength[g] = 0.0
end
--src[1] = 1.0
mg_src = lbs.VolumetricSource.Create({ block_ids = { 0 }, group_strength = strength })

-- Setup Physics
pquad0 = aquad.CreateGLProductQuadrature1DSlab(80)
lbs_block = {
  mesh = grid,
  num_groups = num_groups,
  groupsets = {
    {
      groups_from_to = { 0, 62 },
      angular_quadrature = pquad0,
      angle_aggregation_num_subsets = 1,
      inner_linear_method = "petsc_gmres",
      l_abs_tol = 1.0e-6,
      l_max_its = 300,
      gmres_restart_interval = 100,
    },
    {
      groups_from_to = { 63, num_groups - 1 },
      angular_quadrature = pquad0,
      angle_aggregation_num_subsets = 1,
      inner_linear_method = "petsc_gmres",
      l_abs_tol = 1.0e-6,
      l_max_its = 300,
      gmres_restart_interval = 100,
    },
  },
  xs_map = {
    { block_ids = { 0 }, xs = xs_3_170 },
  },
  sweep_type = "CBC",
}

bsrc = {}
for g = 1, num_groups do
  bsrc[g] = 0.0
end
bsrc[1] = 1.0 / 2

lbs_options = {
  boundary_conditions = {
    {
      name = "zmin",
      type = "isotropic",
      group_strength = bsrc,
    },
  },
  scattering_order = 5,
  save_angular_flux = true,
  max_ags_iterations = 1,
  volumetric_sources = { mg_src },
}

phys1 = lbs.DiscreteOrdinatesSolver.Create(lbs_block)
phys1:SetOptions(lbs_options)

-- Initialize and Execute Solver
ss_solver = lbs.SteadyStateSolver.Create({ lbs_solver = phys1 })

ss_solver:Initialize()
ss_solver:Execute()

-- Get field functions
fflist = lbs.GetScalarFieldFunctionList(phys1)

-- Line plot
--Testing consolidated interpolation
cline = fieldfunc.FieldFunctionInterpolationLine.Create()
cline:SetInitialPoint({ x = 0.0, y = 0.0, z = 0.0001 + xmin })
cline:SetFinalPoint({ x = 0.0, y = 0.0, z = 29.999 + xmin })
cline:SetNumberOfPoints(50)

for k = 165, 165 do
  cline:AddFieldFunction(fflist[k])
end

cline:Initialize()
cline:Execute()

-- Volume integrations
vol0 = logvol.RPPLogicalVolume.Create({ infx = true, infy = true, infz = true })
ffi1 = fieldfunc.FieldFunctionInterpolationVolume.Create()
curffi = ffi1
curffi:SetOperationType(OP_MAX)
curffi:SetLogicalVolume(vol0)
curffi:AddFieldFunction(fflist[1])

curffi:Initialize()
curffi:Execute()
maxval = curffi:GetValue()

log.Log(LOG_0, string.format("Max-value1=%.5f", maxval))

ffi2 = fieldfunc.FieldFunctionInterpolationVolume.Create()
curffi = ffi2
curffi:SetOperationType(OP_MAX)
curffi:SetLogicalVolume(vol0)
curffi:AddFieldFunction(fflist[160])

curffi:Initialize()
curffi:Execute()
maxval = curffi:GetValue()

log.Log(LOG_0, string.format("Max-value2=%.5e", maxval))

-- Exports
if master_export == nil then
  fieldfunc.ExportToCSV(cline)
end

-- Plots
if location_id == 0 and master_export == nil then
  local handle = io.popen("python3 ZLFFI00.py")
end
