--- Final k-eigenvalue    :         1.1925596 (265)

dofile("mesh/gmesh_coarse.lua")
dofile("materials/materials.lua")

-- Setup Physics

-- Angular quadrature
pquad = aquad.CreateProductQuadrature(GAUSS_LEGENDRE_CHEBYSHEV, 2, 2)
aquad.OptimizeForPolarSymmetry(pquad, 4.0 * math.pi)

-- Solver
if string.find(k_method, "scdsa") or string.find(k_method, "smm") then
  inner_linear_method = "richardson"
  l_max_its = 1
else
  inner_linear_method = "gmres"
  l_max_its = 5
end

phys1 = lbs.DiscreteOrdinatesSolver.Create({
  num_groups = num_groups,
  groupsets = {
    {
      groups_from_to = { 0, num_groups - 1 },
      angular_quadrature_handle = pquad,
      inner_linear_method = inner_linear_method,
      l_max_its = l_max_its,
      l_abs_tol = 1.0e-10,
      angle_aggregation_type = "polar",
      angle_aggregation_num_subsets = 1,
      groupset_num_subsets = 1,
    },
  },
  options = {
    boundary_conditions = {
      { name = "xmin", type = "reflecting" },
      { name = "ymin", type = "reflecting" },
    },
    scattering_order = 1,
    verbose_outer_iterations = true,
    verbose_inner_iterations = true,
    power_field_function_on = true,
    power_default_kappa = 1.0,
    power_normalization = 1.0,
    save_angular_flux = true,
  },
  sweep_type = "AAH",
})

-- Execute Solver
if k_method == "pi" then
  k_solver = lbs.PowerIterationKEigen.Create({
    lbs_solver_handle = phys1,
    k_tol = 1.0e-8,
  })
elseif k_method == "pi_scdsa" then
  k_solver = lbs.PowerIterationKEigenSCDSA.Create({
    lbs_solver_handle = phys1,
    diff_accel_sdm = "pwld",
    accel_pi_verbose = true,
    k_tol = 1.0e-8,
    accel_pi_k_tol = 1.0e-8,
    accel_pi_max_its = 50,
  })
elseif k_method == "pi_scdsa_pwlc" then
  k_solver = lbs.PowerIterationKEigenSCDSA.Create({
    lbs_solver_handle = phys1,
    diff_accel_sdm = "pwlc",
    accel_pi_verbose = true,
    k_tol = 1.0e-8,
    accel_pi_k_tol = 1.0e-8,
    accel_pi_max_its = 50,
  })
elseif k_method == "pi_smm" then
  k_solver = lbs.PowerIterationKEigenSMM.Create({
    lbs_solver_handle = phys1,
    accel_pi_verbose = true,
    k_tol = 1.0e-8,
    accel_pi_k_tol = 1.0e-8,
    accel_pi_max_its = 30,
    diff_sdm = "pwlc",
  })
elseif k_method == "pi_smm_pwld" then
  k_solver = lbs.PowerIterationKEigenSMM.Create({
    lbs_solver_handle = phys1,
    accel_pi_verbose = true,
    k_tol = 1.0e-8,
    accel_pi_k_tol = 1.0e-8,
    accel_pi_max_its = 30,
    diff_sdm = "pwld",
  })
elseif k_method == "jfnk" then
  k_solver = lbs.NonLinearKEigen.Create({
    lbs_solver_handle = phys1,
    nl_max_its = 50,
    nl_abs_tol = 1.0e-10,
    nl_rel_tol = 1.0e-10,
    l_max_its = 20,
    num_initial_power_iterations = 2,
  })
else
  log.Log(
    LOG_0ERROR,
    'k_method must be specified. "pi", '
      .. '"pi_scdsa", "pi_scdsa_pwlc", "pi_smm", "pi_smm_pwld", '
      .. 'or "jfnk"'
  )
  return
end

solver.Initialize(k_solver)
solver.Execute(k_solver)

if master_export == nil then
  fflist, count = lbs.GetScalarFieldFunctionList(phys1)
  fieldfunc.ExportToVTKMulti(fflist, "solutions/ZPhi")
end
