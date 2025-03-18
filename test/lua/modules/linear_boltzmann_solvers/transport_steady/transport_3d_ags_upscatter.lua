-- 3D Transport test with multiple groupsets and groupset-to-groupset upscattering
-- SDM: PWLD

num_procs = 4

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
meshgen1 = mesh.ExtruderMeshGenerator.Create({
  inputs = {
    mesh.FromFileMeshGenerator.Create({
      filename = "../../../../assets/mesh/TriangleMesh2x2Cuts.obj",
    }),
  },
  layers = { { z = 0.4, n = 2 }, { z = 0.8, n = 2 }, { z = 1.2, n = 2 }, { z = 1.6, n = 2 } }, -- layers
  partitioner = mesh.KBAGraphPartitioner.Create({
    nx = 2,
    ny = 2,
    xcuts = { 0.0 },
    ycuts = { 0.0 },
  }),
})
grid = meshgen1:Execute()

-- Set block IDs
vol0 = logvol.RPPLogicalVolume.Create({ infx = true, infy = true, infz = true })
grid:SetBlockIDFromLogicalVolume(vol0, 0, true)

num_groups = 3
xs_upscatter = xs.LoadFromOpenSn("simple_upscatter.xs")

strength = {}
for g = 1, num_groups do
  strength[g] = 1.0
end
mg_src = lbs.VolumetricSource.Create({ block_ids = { 0 }, group_strength = strength })

-- Setup Physics
pquad0 = aquad.CreateGLCProductQuadrature3DXYZ(4, 8)

lbs_block = {
  mesh = grid,
  num_groups = num_groups,
  groupsets = {
    {
      groups_from_to = { 0, 0 },
      angular_quadrature = pquad0,
      inner_linear_method = "petsc_gmres",
      l_abs_tol = 1.0e-6,
      l_max_its = 300,
      gmres_restart_interval = 30,
    },
    {
      groups_from_to = { 1, 1 },
      angular_quadrature = pquad0,
      inner_linear_method = "petsc_gmres",
      l_abs_tol = 1.0e-6,
      l_max_its = 300,
      gmres_restart_interval = 30,
    },
    {
      groups_from_to = { 2, 2 },
      angular_quadrature = pquad0,
      inner_linear_method = "petsc_gmres",
      l_abs_tol = 1.0e-6,
      l_max_its = 300,
      gmres_restart_interval = 30,
    },
  },
  xs_map = {
    { block_ids = { 0 }, xs = xs_upscatter },
  },
}

lbs_options = {
  scattering_order = 0,
  verbose_ags_iterations = true,
  max_ags_iterations = 30,
  ags_tolerance = 1.0e-6,
  volumetric_sources = { mg_src },
}

phys1 = lbs.DiscreteOrdinatesSolver.Create(lbs_block)
phys1:SetOptions(lbs_options)

-- Initialize and Execute Solver
ss_solver = lbs.SteadyStateSolver.Create({ lbs_problem = phys1 })
ss_solver:Initialize()
ss_solver:Execute()

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
log.Log(LOG_0, string.format("Max-value1=%.5e", maxval))

ffi1 = fieldfunc.FieldFunctionInterpolationVolume.Create()
curffi = ffi1
curffi:SetOperationType(OP_MAX)
curffi:SetLogicalVolume(vol0)
curffi:AddFieldFunction(fflist[2])
curffi:Initialize()
curffi:Execute()
maxval = curffi:GetValue()
log.Log(LOG_0, string.format("Max-value2=%.5e", maxval))

ffi1 = fieldfunc.FieldFunctionInterpolationVolume.Create()
curffi = ffi1
curffi:SetOperationType(OP_MAX)
curffi:SetLogicalVolume(vol0)
curffi:AddFieldFunction(fflist[3])
curffi:Initialize()
curffi:Execute()
maxval = curffi:GetValue()
log.Log(LOG_0, string.format("Max-value3=%.5e", maxval))
