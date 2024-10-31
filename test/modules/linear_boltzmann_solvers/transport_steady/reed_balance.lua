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
mesh.MeshGenerator.Execute(meshgen)

-- Set Material IDs
z_min = 0.0
z_max = widths[1]
for imat = 1, Nmat do
  z_max = z_min + widths[imat]
  log.Log(LOG_0, "imat=" .. imat .. ", zmin=" .. z_min .. ", zmax=" .. z_max)
  lv = logvol.RPPLogicalVolume.Create({ infx = true, infy = true, zmin = z_min, zmax = z_max })
  mesh.SetMaterialIDFromLogicalVolume(lv, imat - 1)
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
  mat.SetProperty(materials[imat], TRANSPORT_XSECTIONS, SIMPLE_ONE_GROUP, total[imat], c[imat])
end

-- Create sources in 1st and 4th materials
mat.SetProperty(materials[1], ISOTROPIC_MG_SOURCE, FROM_ARRAY, { 50. })
mat.SetProperty(materials[4], ISOTROPIC_MG_SOURCE, FROM_ARRAY, { 1. })

-- Angular Quadrature
gl_quad = aquad.CreateProductQuadrature(GAUSS_LEGENDRE, 64)

-- LBS block option
num_groups = 1
lbs_block = {
  num_groups = num_groups,
  groupsets = {
    {
      groups_from_to = { 0, num_groups - 1 },
      angular_quadrature_handle = gl_quad,
      inner_linear_method = "petsc_gmres",
      l_abs_tol = 1.0e-9,
      l_max_its = 300,
      gmres_restart_interval = 30,
    },
  },
  options = {
    scattering_order = 0,
    spatial_discretization = "pwld",
    boundary_conditions = { { name = "zmin", type = "vacuum" }, { name = "zmax", type = "vacuum" } },
  },
}

phys = lbs.DiscreteOrdinatesSolver.Create(lbs_block)

-- Initialize and execute solver
ss_solver = lbs.SteadyStateSolver.Create({ lbs_solver_handle = phys })

solver.Initialize(ss_solver)
solver.Execute(ss_solver)

-- compute particle balance
lbs.ComputeBalance(phys)
