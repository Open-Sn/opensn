-- 3D 2G Wedge import test. Imports EXODUSII.
-- Uses KEigenvalue::Solver with Power Iteration
-- Test: Final k-eigenvalue: 0.9293376

-- Set and check number of processors
num_procs = 4
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
meshgen1 = mesh.MeshGenerator.Create({
    inputs = {
      mesh.FromFileMeshGenerator.Create({
        filename = "fuel_cube.e",
      }),
    },
    partitioner = mesh.PETScGraphPartitioner.Create({ type = "parmetis" }),
})
mesh.MeshGenerator.Execute(meshgen1)

-- Set Materials (Fuel)
materials = {}
materials[1] = mat.AddMaterial("Fuel")
mat.SetProperty(materials[1], TRANSPORT_XSECTIONS, OPENSN_XSFILE, "xs_fuel_g2.xs")

-- Set Quadrature
nazimu = 4
npolar = 4
pquad = aquad.CreateProductQuadrature(GAUSS_LEGENDRE_CHEBYSHEV, nazimu, npolar)

num_groups = 2
-- Initialize the LBSolver
lbs_block = {
    num_groups = num_groups,
    groupsets = {
        {
            groups_from_to = { 0, 1 },
            angular_quadrature_handle = pquad,
            angle_aggregation_type = 'single',
            groupset_num_subsets = 1,
            inner_linear_method = "petsc_gmres",
            l_abs_tol = 1.0e-6,
            l_max_its = 300,
            gmres_restart_interval = 30
        },
    },
}

lbs_options = {
    -- boundary_conditions = {
    --     {name = "xmin", type = "reflecting"},
    --     {name = "xmax", type = "reflecting"},
    --     {name = "ymin", type = "reflecting"},
    --     {name = "ymax", type = "reflecting"},
    --     {name = "zmin", type = "reflecting"},
    --     {name = "zmax", type = "reflecting"},
    -- },
    scattering_order = 1,
    use_precursors = false,
    verbose_inner_iterations = false,
    verbose_outer_iterations = true, 
}
  
phys1 = lbs.DiscreteOrdinatesSolver.Create(lbs_block)
lbs.SetOptions(phys1, lbs_options)

k_solver0 = lbs.PowerIterationKEigen.Create({
     lbs_solver_handle = phys1,
     k_tol = 1e-6,
    })
solver.Initialize(k_solver0)
solver.Execute(k_solver0)

fflist, count = lbs.GetScalarFieldFunctionList(phys1)
fieldfunc.ExportToVTKMulti(fflist, 'flux')
mesh.ExportToPVTU("mesh")
