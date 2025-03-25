# 1D Transport leakage test
# Unit angular flux left boundary condition in a pure absorber with unit
# length and a unit absorption cross section. The analytic solution is:
# j^+ = \int_{0}^{1} \mu e^{-1/\mu} d\mu = 0.10969

# Check num_procs
num_procs = 3
if check_num_procs == None and number_of_processes ~= num_procs then
  log.Log(
    LOG_0ERROR,
    "Incorrect amount of processors. "
      + "Expected "
      + tostring(num_procs)
      + ". Pass check_num_procs=False to override if possible."
  )
  os.exit(False)
end

# Setup mesh
N = 100
L = 1.0
nodes = {}
for i = 1, (N + 1) do
  k = i - 1
  nodes[i] = (i - 1) * L / N
end

meshgen = mesh.OrthogonalMeshGenerator.Create({ node_sets = { nodes } })
grid = meshgen:Execute()
grid:SetUniformBlockID(0)

# Add materials
num_groups = 1
sigma_t = 1.0

xs1g = xs.CreateSimpleOneGroup(sigma_t, 0.0)

# Setup Physics
pquad = aquad.CreateGLProductQuadrature1DSlab(256)
lbs_block = {
  mesh = grid,
  num_groups = num_groups,
  groupsets = {
    {
      groups_from_to = { 0, num_groups - 1 },
      angular_quadrature = pquad,
      angle_aggregation_num_subsets = 1,
      inner_linear_method = "petsc_gmres",
      l_abs_tol = 1.0e-6,
      l_max_its = 300,
      gmres_restart_interval = 100,
    },
  },
  xs_map = {
    { block_ids = { 0 }, xs = xs1g },
  },
}

bsrc = {}
for g = 1, num_groups do
  bsrc[g] = 0.0
end
bsrc[1] = 1.0

lbs_options = {
  boundary_conditions = {
    {
      name = "zmin",
      type = "isotropic",
      group_strength = bsrc,
    },
  },
  scattering_order = 0,
  save_angular_flux = True,
}

phys = lbs.DiscreteOrdinatesSolver.Create(lbs_block)
phys:SetOptions(lbs_options)

ss_solver = lbs.SteadyStateSolver.Create({ lbs_solver = phys })

# Solve the problem
ss_solver:Initialize()
ss_solver:Execute()

# Compute the leakage
leakage = lbs.ComputeLeakage(phys, {})
for k, v in pairs(leakage) do
  log.Log(LOG_0, string.format("%s=%.5e", k, v[1]))
end
