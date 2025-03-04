-- SDM: PWLD

Ng = 64

Npolar = 4
Nazimuthal = 16

meshgen1 = mesh.MeshGenerator.Create({
  inputs = {
    mesh.FromFileMeshGenerator.Create({
      filename = "../../../assets/mesh/InclusionsGmshV2.msh",
    }),
  },
})
grid = meshgen1:Execute()

-- Material
vol0 = logvol.RPPLogicalVolume.Create({ infx = true, infy = true, infz = true })
xs_diag = xs.LoadFromOpenSn("diag_XS_64g_1mom_c0.99.xs")
strength = {}
for g = 1, Ng do
  strength[g] = 0.0
end
strength[1] = 100.0
mg_src = lbs.VolumetricSource.Create({ block_ids = { 1 }, group_strength = strength })

lbs_options = {
  boundary_conditions = {
    { name = "xmin", type = "reflecting" },
    { name = "ymin", type = "reflecting" },
  },
  scattering_order = 0,
  volumetric_sources = { mg_src },
}

-- Quadrature
pquad0 = aquad.CreateGLCProductQuadrature3DXYZ(Npolar, Nazimuthal)

-- Set up solver
gs1 = { 0, Ng - 1 }
lbs_block = {
  mesh = grid,
  num_groups = Ng,
  groupsets = {
    {
      groups_from_to = gs1,
      angular_quadrature = pquad0,
      angle_aggregation_type = "single",
      angle_aggregation_num_subsets = 1,
      inner_linear_method = "petsc_gmres",
      l_abs_tol = 1.0e-6,
      l_max_its = 100,
    },
  },
  xs_map = {
    { block_ids = { 0, 1 }, xs = xs_diag },
  },
}
phys = lbs.DiscreteOrdinatesSolver.Create(lbs_block)
phys:SetOptions(lbs_options)
ss_solver = lbs.SteadyStateSolver.Create({ lbs_solver = phys })

-- Solve
ss_solver:Initialize()
ss_solver:Execute()

fflist = lbs.GetScalarFieldFunctionList(phys)
ffi1 = fieldfunc.FieldFunctionInterpolationVolume.Create()
curffi = ffi1
curffi:SetOperationType(OP_MAX)
curffi:SetLogicalVolume(vol0)
curffi:AddFieldFunction(fflist[1])
curffi:Initialize()
curffi:Execute()
maxval = curffi:GetValue()
log.Log(LOG_0, string.format("Max-value1=%.5f", maxval))
