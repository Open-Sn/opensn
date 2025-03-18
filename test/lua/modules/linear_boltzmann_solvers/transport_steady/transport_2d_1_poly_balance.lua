-- 2D Transport test. Pure scatterer. Balance
num_procs = 1

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
N = 20
L = 5
xmin = -L / 2
dx = L / N
for i = 1, (N + 1) do
  k = i - 1
  nodes[i] = xmin + k * dx
end

meshgen1 = mesh.OrthogonalMeshGenerator.Create({ node_sets = { nodes, nodes } })
grid = meshgen1:Execute()

-- Set block IDs
vol0 = logvol.RPPLogicalVolume.Create({ infx = true, infy = true, infz = true })
grid:SetBlockIDFromLogicalVolume(vol0, 0, true)
vol1 = logvol.RPPLogicalVolume.Create({ xmin = -1000.0, xmax = 0.0, infy = true, infz = true })
grid:SetBlockIDFromLogicalVolume(vol1, 1, true)

num_groups = 1
xs_1g = xs.CreateSimpleOneGroup(1.0, 1.0)

strength = {}
for g = 1, num_groups do
  strength[g] = 0.0
end
mg_src0 = lbs.VolumetricSource.Create({ block_ids = { 0 }, group_strength = strength })
strength[1] = 1.0
mg_src1 = lbs.VolumetricSource.Create({ block_ids = { 1 }, group_strength = strength })

-- Setup Physics
fac = 1
pquad = aquad.CreateGLCProductQuadrature2DXY(6 * fac, 16 * fac)

lbs_block = {
  mesh = grid,
  num_groups = num_groups,
  groupsets = {
    {
      groups_from_to = { 0, 0 },
      angular_quadrature = pquad,
      angle_aggregation_num_subsets = 1,
      inner_linear_method = "petsc_gmres",
      l_abs_tol = 1.0e-8,
      l_max_its = 300,
      gmres_restart_interval = 100,
    },
  },
  xs_map = {
    { block_ids = { 0, 1 }, xs = xs_1g },
  },
}

lbs_options = {
  scattering_order = 0,
  volumetric_sources = { mg_src0, mg_src1 },
}

phys1 = lbs.DiscreteOrdinatesProblem.Create(lbs_block)
phys1:SetOptions(lbs_options)

ss_solver = lbs.SteadyStateSolver.Create({ lbs_problem = phys1 })

ss_solver:Initialize()
ss_solver:Execute()

phys1:ComputeBalance()

-- Get field functions
fflist = lbs.GetScalarFieldFunctionList(phys1)

-- Volume integrations
ffi1 = fieldfunc.FieldFunctionInterpolationVolume.Create()
curffi = ffi1
curffi:SetOperationType(OP_MAX)
curffi:SetLogicalVolume(vol0)
curffi:AddFieldFunction(fflist[1])

curffi:Initialize()
curffi:Execute()
maxval = curffi:GetValue()

log.Log(LOG_0, string.format("Max-value1=%.5f", maxval))

-- Volume integrations
ffi1 = fieldfunc.FieldFunctionInterpolationVolume.Create()
curffi = ffi1
curffi:SetOperationType(OP_MAX)
curffi:SetLogicalVolume(vol0)
curffi:AddFieldFunction(fflist[1])

curffi:Initialize()
curffi:Execute()
maxval = curffi:GetValue()

log.Log(LOG_0, string.format("Max-value2=%.5e", maxval))

-- Exports
if master_export == nil then
  fieldfunc.ExportToVTK(fflist[1], "ZPhi3D", "Phi")
end
