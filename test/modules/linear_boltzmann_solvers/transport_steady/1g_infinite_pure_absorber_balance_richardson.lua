-- Infinite, 1-group, pure absorber with balance
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

-- Set Material IDs
grid:SetUniformMaterialID(0)

materials = {}
materials[1] = mat.AddMaterial("TestMat")

num_groups = 1

-- Add cross sections to materials
xs1g = xs.CreateSimpleOneGroup(1.0, 0.0)
materials[1]:SetTransportXSections(xs1g)

src = {}
src[1] = 1.0
mg_src = xs.IsotropicMultiGroupSource.FromArray(src)
materials[1]:SetIsotropicMGSource(mg_src)

-- Angular Quadrature
pquad = aquad.CreateProductQuadrature(GAUSS_LEGENDRE_CHEBYSHEV, 2, 2)

-- LBS block option
lbs_block = {
  mesh = grid,
  num_groups = num_groups,
  groupsets = {
    {
      groups_from_to = { 0, num_groups - 1 },
      angular_quadrature = pquad,
      inner_linear_method = "petsc_richardson",
      l_abs_tol = 1.0e-9,
      l_max_its = 300,
    },
  },
  options = {
    boundary_conditions = {
      { name = "xmin", type = "reflecting" },
      { name = "xmax", type = "reflecting" },
      { name = "ymin", type = "reflecting" },
      { name = "ymax", type = "reflecting" },
      { name = "zmin", type = "reflecting" },
      { name = "zmax", type = "reflecting" },
    },
  },
}

phys = lbs.DiscreteOrdinatesSolver.Create(lbs_block)

-- Initialize and execute solver
ss_solver = lbs.SteadyStateSolver.Create({ lbs_solver = phys })

ss_solver:Initialize()
ss_solver:Execute()

-- compute particle balance
phys:ComputeBalance()
