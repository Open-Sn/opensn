-- 3D 172G KEigenvalue::Solver test using power iteration and OpenMC MGXS library
-- Test: Final k-eigenvalue: 1.5029618
num_procs = 4

-- Cells
Nx = 5
Ny = 5
Nz = 5

-- Dimensions
Lx = 2.0
Ly = 2.0
Lz = 2.0

-- Mesh
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

-- Materials

materials = {}
materials[1] = mat.AddMaterial("Fissile Material")
xs_uo2 = xs.LoadFromOpenMC("uo2.h5", "set1", 294.0)
materials[1]:SetTransportXSections(xs_uo2)
grid:SetUniformMaterialID(0)

-- Solver

num_groups = 172
lbs_block = {
  mesh = grid,
  num_groups = num_groups,
  groupsets = {
    {
      groups_from_to = { 0, num_groups - 1 },
      angular_quadrature = aquad.CreateGLCProductQuadrature3DXYZ(4, 32),
      inner_linear_method = "classic_richardson",
      l_max_its = 500,
      l_abs_tol = 1.0e-12,
    },
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

k_solver0 = lbs.NonLinearKEigen.Create({
  lbs_solver = phys,
  nl_max_its = 500,
  nl_abs_tol = 1.0e-8,
})
k_solver0:Initialize()
k_solver0:Execute()
