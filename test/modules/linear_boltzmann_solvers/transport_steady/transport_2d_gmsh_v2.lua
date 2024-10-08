-- SDM: PWLD

Ng = 64

Npolar = 4
Nazimuthal = 2

meshgen1 = mesh.MeshGenerator.Create({
  inputs = {
    mesh.FromFileMeshGenerator.Create({
      filename = "../../../assets/mesh/Rectangular2D2MatGmshV2.msh",
    }),
  },
})
mesh.MeshGenerator.Execute(meshgen1)

-- Material
vol0 = logvol.RPPLogicalVolume.Create({ infx = true, infy = true, infz = true })
materials = {}
materials[0] = mat.AddMaterial("Test Material")
materials[1] = mat.AddMaterial("Test Material")
mat.SetProperty(materials[0], TRANSPORT_XSECTIONS, OPENSN_XSFILE, "diag_XS_64g_1mom_c0.99.xs")
mat.SetProperty(materials[1], TRANSPORT_XSECTIONS, OPENSN_XSFILE, "diag_XS_64g_1mom_c0.99.xs")
src = {}
for g = 1, Ng do
  src[g] = 0.0
end
src[1] = 100.0
mat.SetProperty(materials[1], ISOTROPIC_MG_SOURCE, FROM_ARRAY, src)

lbs_options = {
  boundary_conditions = {
    { name = "xmin", type = "reflecting" },
    { name = "ymin", type = "reflecting" },
  },
  scattering_order = 0,
}

-- Quadrature
pquad0 = aquad.CreateProductQuadrature(GAUSS_LEGENDRE_CHEBYSHEV, Npolar, Nazimuthal)

-- Set up solver
gs1 = { 0, Ng - 1 }
lbs_block = {
  num_groups = Ng,
  groupsets = {
    {
      groups_from_to = gs1,
      angular_quadrature_handle = pquad0,
      angle_aggregation_type = "single",
      angle_aggregation_num_subsets = 1,
      groupset_num_subsets = 1,
      inner_linear_method = "gmres",
      l_abs_tol = 1.0e-6,
      l_max_its = 300,
    },
  },
}
phys = lbs.DiscreteOrdinatesSolver.Create(lbs_block)
lbs.SetOptions(phys, lbs_options)
ss_solver = lbs.SteadyStateSolver.Create({ lbs_solver_handle = phys })

-- Solve
solver.Initialize(ss_solver)
solver.Execute(ss_solver)

fflist, count = lbs.GetScalarFieldFunctionList(phys)
ffi1 = fieldfunc.FFInterpolationCreate(VOLUME)
curffi = ffi1
fieldfunc.SetProperty(curffi, OPERATION, OP_MAX)
fieldfunc.SetProperty(curffi, LOGICAL_VOLUME, vol0)
fieldfunc.SetProperty(curffi, ADD_FIELDFUNCTION, fflist[1])
fieldfunc.Initialize(curffi)
fieldfunc.Execute(curffi)
maxval = fieldfunc.GetValue(curffi)
log.Log(LOG_0, string.format("Max-value1=%.5f", maxval))
