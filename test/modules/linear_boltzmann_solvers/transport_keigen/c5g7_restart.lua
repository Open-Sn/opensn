--- Final k-eigenvalue    :         1.1925596 (265)

-- Mesh
meshgen1 = mesh.MeshGenerator.Create({
  inputs = {
    mesh.FromFileMeshGenerator.Create({
      filename = "c5g7/mesh/2D_c5g7_coarse.msh",
    }),
  },
})
mesh.MeshGenerator.Execute(meshgen1)

-- Materials
xss = {}

for m = 0, 6 do
  xss[tostring(m)] = xs.Create()
end

xs.Set(xss["0"], OPENSN_XSFILE, "c5g7/materials/XS_water.xs")
xs.Set(xss["1"], OPENSN_XSFILE, "c5g7/materials/XS_UO2.xs")
xs.Set(xss["2"], OPENSN_XSFILE, "c5g7/materials/XS_7pMOX.xs")
xs.Set(xss["3"], OPENSN_XSFILE, "c5g7/materials/XS_guide_tube.xs")
xs.Set(xss["4"], OPENSN_XSFILE, "c5g7/materials/XS_4_3pMOX.xs")
xs.Set(xss["5"], OPENSN_XSFILE, "c5g7/materials/XS_8_7pMOX.xs")
xs.Set(xss["6"], OPENSN_XSFILE, "c5g7/materials/XS_fission_chamber.xs")
water_xs = xs.Get(xss["0"])
num_groups = water_xs["num_groups"]
log.Log(LOG_0, "Num groups: " .. tostring(num_groups))

materials = {}
for m = 0, 6 do
  key = tostring(m)
  materials[key] = mat.AddMaterial("Material_" .. key)
  mat.SetProperty(materials[key], TRANSPORT_XSECTIONS, EXISTING, xss[key])
end

-- Angular quadrature
pquad = aquad.CreateProductQuadrature(GAUSS_LEGENDRE_CHEBYSHEV, 2, 2)
aquad.OptimizeForPolarSymmetry(pquad, 4.0 * math.pi)

-- Solver
phys1 = lbs.DiscreteOrdinatesSolver.Create({
  num_groups = num_groups,
  groupsets = {
    {
      groups_from_to = { 0, num_groups - 1 },
      angular_quadrature_handle = pquad,
      inner_linear_method = "gmres",
      l_max_its = 5,
      l_abs_tol = 1.0e-10,
      angle_aggregation_type = "polar",
      angle_aggregation_num_subsets = 1,
      groupset_num_subsets = 1,
    },
  },
  options = {
    boundary_conditions = {
      { name = "xmax", type = "reflecting" },
      { name = "ymax", type = "reflecting" },
    },
    scattering_order = 1,
    verbose_outer_iterations = true,
    verbose_inner_iterations = true,
    power_field_function_on = true,
    power_default_kappa = 1.0,
    power_normalization = 1.0,
    save_angular_flux = true,
    read_restart_path = "c5g7_restart/c5g7",
  },
  sweep_type = "CBC",
})

k_solver = lbs.PowerIterationKEigen.Create({
  lbs_solver_handle = phys1,
  k_tol = 1.0e-8,
})

solver.Initialize(k_solver)
solver.Execute(k_solver)

if master_export == nil then
  fflist, count = lbs.GetScalarFieldFunctionList(phys1)
  fieldfunc.ExportToVTKMulti(fflist, "solutions/ZPhi")
end
