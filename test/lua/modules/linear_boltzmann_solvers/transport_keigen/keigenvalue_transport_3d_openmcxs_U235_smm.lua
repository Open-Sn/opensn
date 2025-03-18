-- 3D 84G keigenvalue test using OpenMC MGXS cross-sections and power iteration with SMM
--
-- Test: Final k-eigenvalue: 2.280431
num_procs = 4

--
-- Mesh
--

-- Cells
Nx = 5
Ny = 5
Nz = 5

-- Dimensions
Lx = 2.0
Ly = 2.0
Lz = 2.0

xmesh = {}
xmin = 0.0
dx = Lx / Nx
for i = 1, (Nx + 1) do
  k = i - 1
  xmesh[i] = xmin + k * dx
end

ymesh = {}
ymin = 0.0
dy = Ly / Ny
for i = 1, (Ny + 1) do
  k = i - 1
  ymesh[i] = ymin + k * dy
end

zmesh = {}
zmin = 0.0
dz = Lz / Nz
for i = 1, (Nz + 1) do
  k = i - 1
  zmesh[i] = zmin + k * dz
end

meshgen1 = mesh.OrthogonalMeshGenerator.Create({ node_sets = { xmesh, ymesh, zmesh } })
grid = meshgen1:Execute()

grid:SetUniformBlockID(0)

xs_u235 = xs.LoadFromOpenMC("u235_84g.h5", "set1", 294.0)

--
-- Solver
--

num_groups = 84
lbs_block = {
  mesh = grid,
  num_groups = num_groups,
  groupsets = {
    {
      groups_from_to = { 0, num_groups - 1 },
      angular_quadrature = aquad.CreateGLCProductQuadrature3DXYZ(2, 4),
      inner_linear_method = "classic_richardson",
      l_max_its = 2,
      l_abs_tol = 1.0e-12,
    },
  },
  xs_map = {
    { block_ids = { 0 }, xs = xs_u235 },
  },
}

lbs_options = {
  boundary_conditions = {
    { name = "xmin", type = "reflecting" },
    { name = "xmax", type = "reflecting" },
    { name = "ymin", type = "reflecting" },
    { name = "ymax", type = "reflecting" },
    { name = "zmin", type = "reflecting" },
    { name = "zmax", type = "reflecting" },
  },
  scattering_order = 1,
  use_precursors = false,
  verbose_inner_iterations = false,
  verbose_outer_iterations = true,
}

phys = lbs.DiscreteOrdinatesSolver.Create(lbs_block)
phys:SetOptions(lbs_options)

k_solver0 = lbs.PowerIterationKEigenSMM.Create({
  lbs_problem = phys,
  accel_pi_verbose = true,
  k_tol = 1.0e-8,
  accel_pi_k_tol = 1.0e-8,
  accel_pi_max_its = 30,
  diff_sdm = "pwld",
})
k_solver0:Initialize()
k_solver0:Execute()
