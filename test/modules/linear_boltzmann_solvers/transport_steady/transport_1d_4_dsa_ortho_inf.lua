-- 1D LinearBSolver test of a block of graphite mimicking an infinite medium. DSA and TG
-- SDM: PWLD
-- Test: WGS groups [0-62] Iteration    22 Residual 9.6079e-08 CONVERGED
-- and   WGS groups [63-167] Iteration    59 Residual 4.73732e-07 CONVERGED
num_procs = 1

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
nodes = {}
N = 1
L = 1e6
--N=10
--L=200e6
xmin = -L / 2
--xmin = 0.0
dx = L / N
for i = 1, (N + 1) do
  k = i - 1
  nodes[i] = xmin + k * dx
end

meshgen1 = mesh.OrthogonalMeshGenerator.Create({ node_sets = { nodes } })
grid = meshgen1:Execute()

-- Set Material IDs
grid:SetUniformMaterialID(0)

-- Add materials
materials = {}
materials[1] = mat.AddMaterial("Test Material")
materials[2] = mat.AddMaterial("Test Material2")

num_groups = 168
xs_graphite = xs.LoadFromOpenSn("xs_graphite_pure.xs")
materials[1]:SetTransportXSections(xs_graphite)
xs_air = xs.LoadFromOpenSn("xs_air50RH.xs")
materials[2]:SetTransportXSections(xs_air)

strength = {}
for g = 1, num_groups do
  strength[g] = 0.0
end
strength[1] = 1.0
mg_src = lbs.VolumetricSource.Create({ block_ids = { 0 }, group_strength = strength })

-- Setup Physics
pquad = aquad.GLProductQuadrature1DSlab.Create({ Npolar = 4 })

lbs_block = {
  mesh = grid,
  num_groups = num_groups,
  groupsets = {
    {
      groups_from_to = { 0, 62 },
      angular_quadrature = pquad,
      angle_aggregation_num_subsets = 1,
      inner_linear_method = "petsc_gmres",
      l_abs_tol = 1.0e-6,
      l_max_its = 1000,
      gmres_restart_interval = 30,
      apply_wgdsa = true,
      wgdsa_l_abs_tol = 1.0e-2,
    },
    {
      groups_from_to = { 63, num_groups - 1 },
      angular_quadrature = pquad,
      angle_aggregation_num_subsets = 1,
      inner_linear_method = "petsc_gmres",
      l_abs_tol = 1.0e-6,
      l_max_its = 1000,
      gmres_restart_interval = 30,
      apply_wgdsa = true,
      apply_tgdsa = true,
      wgdsa_l_abs_tol = 1.0e-2,
    },
  },
}

lbs_options = {
  scattering_order = 1,
  volumetric_sources = { mg_src },
}

phys1 = lbs.DiscreteOrdinatesSolver.Create(lbs_block)
phys1:SetOptions(lbs_options)

-- Initialize and Execute Solver
ss_solver = lbs.SteadyStateSolver.Create({ lbs_solver = phys1 })

ss_solver:Initialize()
ss_solver:Execute()

-- Get field functions
fflist = lbs.GetScalarFieldFunctionList(phys1)

-- Exports
if master_export == nil then
  fieldfunc.ExportToVTKMulti(fflist, "ZPhi")
end

-- Plots
