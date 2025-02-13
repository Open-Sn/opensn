-- 1D 1G KEigenvalue::Solver test using power iteration
-- Test: Final k-eigenvalue: 0.9995433
num_procs = 4

-- NOTE: For command line inputs, specify as:
--       variable=[[argument]]

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

MPIBarrier()

-- ##################################################
-- ##### Parameters #####
-- ##################################################

-- Mesh variables
if L == nil then
  L = 100.0
end
if n_cells == nil then
  n_cells = 50
end

-- Transport angle information
if n_angles == nil then
  n_angles = 16
end
if scat_order == nil then
  scat_order = 0
end

-- k-eigenvalue iteration parameters
if kes_max_iterations == nil then
  kes_max_iterations = 5000
end
if kes_tolerance == nil then
  kes_tolerance = 1e-8
end

-- Source iteration parameters
if si_max_iterations == nil then
  si_max_iterations = 500
end
if si_tolerance == nil then
  si_tolerance = 1e-8
end

-- Delayed neutrons
if use_precursors == nil then
  use_precursors = true
end

-- ##################################################
-- ##### Run problem #####
-- ##################################################

-- Setup mesh
nodes = {}
dx = L / n_cells
for i = 0, n_cells do
  nodes[i + 1] = i * dx
end

meshgen1 = mesh.OrthogonalMeshGenerator.Create({ node_sets = { nodes } })
grid = meshgen1:Execute()

-- Set Material IDs
grid:SetUniformMaterialID(0)

-- Add materials
materials = {}
materials[1] = mat.AddMaterial("Fissile Material")

xs_simple_fissile = xs.LoadFromOpenSn("simple_fissile.xs")
materials[1]:SetTransportXSections(xs_simple_fissile)

-- Setup Physics
num_groups = 1
lbs_block = {
  mesh = grid,
  num_groups = num_groups,
  groupsets = {
    {
      groups_from_to = { 0, num_groups - 1 },
      angular_quadrature = aquad.CreateProductQuadrature(GAUSS_LEGENDRE, n_angles, -1),
      inner_linear_method = "petsc_gmres",
      l_max_its = si_max_iterations,
      l_abs_tol = si_tolerance,
    },
  },
}

lbs_options = {
  scattering_order = scat_order,

  use_precursors = use_precursors,

  verbose_inner_iterations = false,
  verbose_outer_iterations = true,
}

phys = lbs.DiscreteOrdinatesSolver.Create(lbs_block)
phys:SetOptions(lbs_options)

k_solver0 = lbs.NonLinearKEigen.Create({
  lbs_solver = phys,
  nl_max_its = kes_max_iterations,
  nl_abs_tol = kes_tolerance,
})
k_solver0:Initialize()
k_solver0:Execute()

-- Get field functions
-- Line plot
-- Volume integrations
-- Exports
-- Plots
