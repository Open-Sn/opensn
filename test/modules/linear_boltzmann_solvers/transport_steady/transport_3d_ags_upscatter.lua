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
      filename = "../../../assets/mesh/TriangleMesh2x2Cuts.obj",
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
mesh.MeshGenerator.Execute(meshgen1)

-- Set Material IDs
vol0 = logvol.RPPLogicalVolume.Create({ infx = true, infy = true, infz = true })
mesh.SetMaterialIDFromLogicalVolume(vol0, 0)

-- Add materials
materials = {}
materials[1] = mat.AddMaterial("Test Material")

num_groups = 3
mat.SetProperty(materials[1], TRANSPORT_XSECTIONS, OPENSN_XSFILE, "simple_upscatter.xs")

src = {}
for g = 1, num_groups do
  src[g] = 1.0
end
mat.SetProperty(materials[1], ISOTROPIC_MG_SOURCE, FROM_ARRAY, src)

-- Setup Physics
pquad0 = aquad.CreateProductQuadrature(GAUSS_LEGENDRE_CHEBYSHEV, 2, 2)

lbs_block = {
  num_groups = num_groups,
  groupsets = {
    {
      groups_from_to = { 0, 0 },
      angular_quadrature_handle = pquad0,
      inner_linear_method = "petsc_gmres",
      l_abs_tol = 1.0e-6,
      l_max_its = 300,
      gmres_restart_interval = 30,
    },
    {
      groups_from_to = { 1, 1 },
      angular_quadrature_handle = pquad0,
      inner_linear_method = "petsc_gmres",
      l_abs_tol = 1.0e-6,
      l_max_its = 300,
      gmres_restart_interval = 30,
    },
    {
      groups_from_to = { 2, 2 },
      angular_quadrature_handle = pquad0,
      inner_linear_method = "petsc_gmres",
      l_abs_tol = 1.0e-6,
      l_max_its = 300,
      gmres_restart_interval = 30,
    },
  },
}

lbs_options = {
  scattering_order = 0,
  verbose_ags_iterations = true,
  max_ags_iterations = 30,
  ags_tolerance = 1.0e-6,
}

phys1 = lbs.DiscreteOrdinatesSolver.Create(lbs_block)
lbs.SetOptions(phys1, lbs_options)

-- Initialize and Execute Solver
ss_solver = lbs.SteadyStateSolver.Create({ lbs_solver_handle = phys1 })
solver.Initialize(ss_solver)
solver.Execute(ss_solver)

-- Get field functions
fflist, count = lbs.GetScalarFieldFunctionList(phys1)

-- Volume integrations
ffi1 = fieldfunc.FFInterpolationCreate(VOLUME)
curffi = ffi1
fieldfunc.SetProperty(curffi, OPERATION, OP_MAX)
fieldfunc.SetProperty(curffi, LOGICAL_VOLUME, vol0)
fieldfunc.SetProperty(curffi, ADD_FIELDFUNCTION, fflist[1])
fieldfunc.Initialize(curffi)
fieldfunc.Execute(curffi)
maxval = fieldfunc.GetValue(curffi)
log.Log(LOG_0, string.format("Max-value1=%.5e", maxval))

ffi1 = fieldfunc.FFInterpolationCreate(VOLUME)
curffi = ffi1
fieldfunc.SetProperty(curffi, OPERATION, OP_MAX)
fieldfunc.SetProperty(curffi, LOGICAL_VOLUME, vol0)
fieldfunc.SetProperty(curffi, ADD_FIELDFUNCTION, fflist[2])
fieldfunc.Initialize(curffi)
fieldfunc.Execute(curffi)
maxval = fieldfunc.GetValue(curffi)
log.Log(LOG_0, string.format("Max-value2=%.5e", maxval))

ffi1 = fieldfunc.FFInterpolationCreate(VOLUME)
curffi = ffi1
fieldfunc.SetProperty(curffi, OPERATION, OP_MAX)
fieldfunc.SetProperty(curffi, LOGICAL_VOLUME, vol0)
fieldfunc.SetProperty(curffi, ADD_FIELDFUNCTION, fflist[3])
fieldfunc.Initialize(curffi)
fieldfunc.Execute(curffi)
maxval = fieldfunc.GetValue(curffi)
log.Log(LOG_0, string.format("Max-value3=%.5e", maxval))
