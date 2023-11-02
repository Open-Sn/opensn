-- 2D 2G KEigenvalue::Solver test using Power Iteration
-- Test: Final k-eigenvalue: 0.5969127

dofile("utils/QBlock_mesh.lua")
dofile("utils/QBlock_materials.lua") --num_groups assigned here

--############################################### Setup Physics
pquad = chiCreateProductQuadrature(GAUSS_LEGENDRE_CHEBYSHEV,4, 4)
chiOptimizeAngularQuadratureForPolarSymmetry(pquad, 4.0*math.pi)

lbs_block =
{
  num_groups = num_groups,
  groupsets =
  {
    {
      groups_from_to = {0, num_groups-1},
      angular_quadrature_handle = pquad,
      inner_linear_method = "richardson",
      l_max_its = 1,
      gmres_restart_interval = 50,
      l_abs_tol = 1.0e-10,
      groupset_num_subsets = 1,
    }
  },
  options =
  {
    boundary_conditions = { { name = "xmin", type = "reflecting"},
                            { name = "ymin", type = "reflecting"} },
    scattering_order = 2,

    use_precursors = false,

    verbose_inner_iterations = true,
    verbose_outer_iterations = true,
    save_angular_flux = true
  }
}

--lbs_options =
--{
--  boundary_conditions = { { name = "xmin", type = "reflecting"},
--                          { name = "ymin", type = "reflecting"} },
--  scattering_order = 2,
--
--  use_precursors = false,
--
--  verbose_inner_iterations = false,
--  verbose_outer_iterations = true,
--}

phys1 = lbs.DiscreteOrdinatesSolver.Create(lbs_block)
--lbs.SetOptions(phys1, lbs_options)


--k_solver0 = lbs.XXPowerIterationKEigen.Create({ lbs_solver_handle = phys1, })
k_solver0 = lbs.XXPowerIterationKEigenSCDSA.Create
({
  lbs_solver_handle = phys1,
  diff_accel_sdm = "pwld",
  accel_pi_verbose = false,
  k_tol = 1.0e-8
})
chiSolverInitialize(k_solver0)
chiSolverExecute(k_solver0)


fflist,count = chiLBSGetScalarFieldFunctionList(phys1)

--chiExportMultiFieldFunctionToVTK(fflist,"tests/BigTests/QBlock/solutions/Flux")

-- Reference value k_eff = 0.5969127