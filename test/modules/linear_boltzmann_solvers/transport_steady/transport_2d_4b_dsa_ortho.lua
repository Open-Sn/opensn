-- 2D LinearBSolver Same as 4a but with reflective BCs. DSA and TG
-- SDM: PWLD
-- Test: WGS groups [0-62] Iteration    54 Residual 5.00021e-07 CONVERGED
-- and   WGS groups [63-167] Iteration    56 Residual 9.73954e-07 CONVERGED
num_procs = 4

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
N = 20
L = 100
--N=10
--L=200e6
xmin = -L / 2
--xmin = 0.0
dx = L / N
for i = 1, (N + 1) do
  k = i - 1
  nodes[i] = xmin + k * dx
end

meshgen1 = mesh.OrthogonalMeshGenerator.Create({ node_sets = { nodes, nodes } })
grid = meshgen1:Execute()

-- Set block IDs
grid:SetUniformBlockID(0)

vol1 = logvol.RPPLogicalVolume.Create({
  xmin = -10.0,
  xmax = 10.0,
  ymin = -10.0,
  ymax = 10.0,
  infz = true,
})
grid:SetMaterialIDFromLogicalVolume(vol1, 1, true)

num_groups = 168
xs_graphite = xs.LoadFromOpenSn("xs_graphite_pure.xs")
xs_air = xs.LoadFromOpenSn("xs_air50RH.xs")

strength = {}
for g = 1, num_groups do
  strength[g] = 0.0
end
strength[1] = 1.0
mg_src0 = lbs.VolumetricSource.Create({ block_ids = { 0 }, group_strength = strength })
strength[1] = 0.0
mg_src1 = lbs.VolumetricSource.Create({ block_ids = { 1 }, group_strength = strength })

-- Setup Physics
pquad0 = aquad.CreateGLCProductQuadrature2DXY(4, 8)

lbs_block = {
  mesh = grid,
  num_groups = num_groups,
  groupsets = {
    {
      groups_from_to = { 0, 62 },
      angular_quadrature = pquad0,
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
      angular_quadrature = pquad0,
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
  xs_map = {
    { block_ids = { 0 }, xs = xs_graphite },
    { block_ids = { 1 }, xs = xs_air },
  },
}

lbs_options = {
  boundary_conditions = {
    { name = "xmin", type = "reflecting" },
    { name = "ymin", type = "reflecting" },
  },
  scattering_order = 1,
  max_ags_iterations = 1,
  volumetric_sources = { mg_src0, mg_src1 },
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
