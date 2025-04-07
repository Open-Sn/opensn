-- Infinite 172-group problem
-- Create Mesh
nodes = {}
N = 2
L = 10
xmin = -L / 2
dx = L / N
for i = 1, (N + 1) do
  k = i - 1
  nodes[i] = xmin + k * dx
end

meshgen = mesh.OrthogonalMeshGenerator.Create({ node_sets = { nodes, nodes, nodes } })
grid = meshgen:Execute()

-- Set block IDs
grid:SetUniformBlockID(0)

num_groups = 172

-- cross sections
xs_hdpe = xs.LoadFromOpenMC("HDPE.h5", "set1", 294.0)

strength = {}
for g = 1, num_groups do
  strength[g] = 0.0
end
strength[1] = 1.0
mg_src = lbs.VolumetricSource.Create({ block_ids = { 0 }, group_strength = strength })

-- Angular Quadrature
pquad = aquad.CreateGLCProductQuadrature3DXYZ(4, 8)

-- LBS block option
lbs_block = {
  mesh = grid,
  num_groups = num_groups,
  groupsets = {
    {
      groups_from_to = { 0, num_groups - 1 },
      angular_quadrature = pquad,
      inner_linear_method = "petsc_gmres",
      l_abs_tol = 1.0e-9,
      l_max_its = 300,
      gmres_restart_interval = 30,
    },
  },
  xs_map = {
    { block_ids = { 0 }, xs = xs_hdpe },
  },
  sweep_type = "CBC",
  options = {
    scattering_order = 0,
    spatial_discretization = "pwld",
    save_angular_flux = true,
    boundary_conditions = {
      { name = "xmin", type = "reflecting" },
      { name = "xmax", type = "reflecting" },
      { name = "ymin", type = "reflecting" },
      { name = "ymax", type = "reflecting" },
      { name = "zmin", type = "reflecting" },
      { name = "zmax", type = "reflecting" },
    },
    volumetric_sources = { mg_src },
  },
}

phys = lbs.DiscreteOrdinatesSolver.Create(lbs_block)

-- Initialize and execute solver
ss_solver = lbs.SteadyStateSolver.Create({ lbs_solver = phys })

ss_solver:Initialize()
ss_solver:Execute()

-- compute particle balance
phys:ComputeBalance()
