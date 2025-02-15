-- Standard Reed 1D 1-group problem
-- Create Mesh
widths = { 2., 1., 2., 1., 2. }
nrefs = { 200, 200, 200, 200, 200 }

Nmat = #widths

nodes = {}
counter = 1
nodes[counter] = 0.
for imat = 1, Nmat do
  dx = widths[imat] / nrefs[imat]
  for i = 1, nrefs[imat] do
    counter = counter + 1
    nodes[counter] = nodes[counter - 1] + dx
  end
end

meshgen = mesh.OrthogonalMeshGenerator.Create({ node_sets = { nodes } })
grid = meshgen:Execute()

-- Set Material IDs
z_min = 0.0
z_max = widths[1]
for imat = 1, Nmat do
  z_max = z_min + widths[imat]
  log.Log(LOG_0, "imat=" .. imat .. ", zmin=" .. z_min .. ", zmax=" .. z_max)
  lv = logvol.RPPLogicalVolume.Create({ infx = true, infy = true, zmin = z_min, zmax = z_max })
  grid:SetMaterialIDFromLogicalVolume(lv, imat - 1, true)
  z_min = z_max
end

-- Create materials
mat_names = { "AbsoSrc", "Abso", "Void", "ScatSrc", "Scat" }
materials = {}
for imat = 1, Nmat do
  materials[imat] = mat.AddMaterial(mat_names[imat])
end

-- Add cross sections to materials
total = { 50., 5., 0., 1., 1. }
c = { 0., 0., 0., 0.9, 0.9 }
for imat = 1, Nmat do
  xs1g = xs.CreateSimpleOneGroup(total[imat], c[imat])
  materials[imat]:SetTransportXSections(xs1g)
end

-- Create sources in 1st and 4th materials
mg_src0 = lbs.VolumetricSource.Create({ block_ids = { 0 }, group_strength = { 50. } })
mg_src1 = lbs.VolumetricSource.Create({ block_ids = { 3 }, group_strength = { 1. } })

-- Angular Quadrature
gl_quad = aquad.CreateProductQuadrature(GAUSS_LEGENDRE, 64, -1)

-- LBS block option
num_groups = 1
lbs_block = {
  mesh = grid,
  num_groups = num_groups,
  groupsets = {
    {
      groups_from_to = { 0, num_groups - 1 },
      angular_quadrature = gl_quad,
      inner_linear_method = "classic_richardson",
      l_abs_tol = 1.0e-9,
      l_max_its = 1000,
    },
  },
  options = {
    scattering_order = 0,
    spatial_discretization = "pwld",
    boundary_conditions = { { name = "zmin", type = "vacuum" }, { name = "zmax", type = "vacuum" } },
    volumetric_sources = { mg_src0, mg_src1 },
  },
}

phys = lbs.DiscreteOrdinatesSolver.Create(lbs_block)

-- Initialize and execute solver
ss_solver = lbs.SteadyStateSolver.Create({ lbs_solver = phys })

ss_solver:Initialize()
ss_solver:Execute()

-- compute particle balance
phys:ComputeBalance()
