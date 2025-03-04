-- 2D 2G KEigenvalue::Solver test using NonLinearK
-- Test: Final k-eigenvalue: 0.5969127

dofile("utils/qblock_mesh.lua")
dofile("utils/qblock_materials.lua") --num_groups assigned here

-- Setup Physics
pquad = aquad.CreateGLCProductQuadrature2DXY(8, 16)

lbs_block = {
  mesh = grid,
  num_groups = num_groups,
  groupsets = {
    {
      groups_from_to = { 0, num_groups - 1 },
      angular_quadrature = pquad,
      inner_linear_method = "petsc_gmres",
      l_max_its = 50,
      gmres_restart_interval = 50,
      l_abs_tol = 1.0e-10,
    },
  },
  xs_map = {
    { block_ids = { 0 }, xs = xss["0"] },
    { block_ids = { 1 }, xs = xss["1"] },
  },
}

lbs_options = {
  boundary_conditions = {
    { name = "xmin", type = "reflecting" },
    { name = "ymin", type = "reflecting" },
  },
  scattering_order = 2,

  use_precursors = false,

  verbose_inner_iterations = false,
  verbose_outer_iterations = true,
}

phys1 = lbs.DiffusionDFEMSolver.Create(lbs_block)
phys1:SetOptions(lbs_options)

k_solver0 = lbs.NonLinearKEigen.Create({ lbs_solver = phys1 })
k_solver0:Initialize()
k_solver0:Execute()

fflist = lbs.GetScalarFieldFunctionList(phys1)

--fieldfunc.ExportToVTKMulti(fflist,"tests/BigTests/QBlock/solutions/Flux")

-- Reference value k_eff = 0.5969127
