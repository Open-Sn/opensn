-- Final k-eigenvalue    :         1.1925596 (265)
dofile("mesh/gmesh_coarse.lua")

--os.exit()

dofile("materials/materials.lua")

--############################################### Setup Physics
--========== ProdQuad
pquad = aquad.CreateProductQuadrature(GAUSS_LEGENDRE_CHEBYSHEV,2, 2)
--pquad = aquad.CreateProductQuadrature(GAUSS_LEGENDRE_CHEBYSHEV,1, 1)
aquad.OptimizeForPolarSymmetry(pquad, 4.0*math.pi)


phys1 = lbs.DiscreteOrdinatesSolver.Create
({
    num_groups = num_groups,
    groupsets =
    {
        {
            groups_from_to = {0, num_groups-1},
            angular_quadrature_handle = pquad,
            inner_linear_method = "gmres",
            l_max_its = 5,
            l_abs_tol = 1.0e-10,
            angle_aggregation_type = "polar",
            angle_aggregation_num_subsets = 1,
            groupset_num_subsets = 1,
        }
    },
    options =
    {
        boundary_conditions = { { name = "xmax", type = "reflecting"},
                                { name = "ymax", type = "reflecting"} },
        scattering_order = 1,

        verbose_inner_iterations = true,
        verbose_outer_iterations = true,

        power_field_function_on = true,
        power_default_kappa = 1.0,
        -- power_normalization = -1.0, --Disabled
        power_normalization = 1.0, --Disabled
        save_angular_flux = true
    },
    sweep_type = "CBC"
})

--############################################### Initialize and Execute Solver
if (k_method == "pi") then
    k_solver = lbs.XXPowerIterationKEigen.Create
    ({
        lbs_solver_handle = phys1,
        k_tol = 1.0e-8
    })
    solver.Initialize(k_solver)
    solver.Execute(k_solver)
elseif (k_method == "pi_scdsa") then
    k_solver = lbs.XXPowerIterationKEigenSCDSA.Create
    ({
        lbs_solver_handle = phys1,
        diff_accel_sdm = "pwld",
        accel_pi_verbose = true,
        k_tol = 1.0e-8,
        accel_pi_k_tol = 1.0e-8,
        accel_pi_max_its = 50,
    })
    LBSGroupsetSetIterativeMethod(phys1, 0, KRYLOV_RICHARDSON_CYCLES)
    LBSGroupsetSetMaxIterations(phys1, 0, 1)

    solver.Initialize(k_solver)
    solver.Execute(k_solver)
elseif (k_method == "pi_scdsa_pwlc") then
    k_solver = lbs.XXPowerIterationKEigenSCDSA.Create
    ({
        lbs_solver_handle = phys1,
        diff_accel_sdm = "pwlc",
        accel_pi_verbose = true,
        k_tol = 1.0e-8,
        accel_pi_k_tol = 1.0e-8,
        accel_pi_max_its = 50,
    })
    LBSGroupsetSetIterativeMethod(phys1, 0, KRYLOV_RICHARDSON_CYCLES)
    LBSGroupsetSetMaxIterations(phys1, 0, 1)

    solver.Initialize(k_solver)
    solver.Execute(k_solver)
elseif (k_method == "jfnk") then
    k_solver = lbs.XXNonLinearKEigen.Create
    ({
        lbs_solver_handle = phys1,
        nl_max_its = 50,
        nl_abs_tol = 1.0e-10,
        nl_rel_tol = 1.0e-10,
        l_max_its = 20,
        num_free_power_iterations = 2,
    })
    solver.Initialize(k_solver)
    solver.Execute(k_solver)
else
    log.Log(LOG_0ERROR, "k_method must be specified. \"pi\", "..
      "\"pi_scdsa\", \"pi_scdsa_pwlc\" or \"jfnk\"");
    os.exit(1)
end

if (master_export == nil) then
    fflist,count = lbs.GetScalarFieldFunctionList(phys1)
    fieldfunc.ExportToVTKMulti(fflist,"solutions/ZPhi")
end
