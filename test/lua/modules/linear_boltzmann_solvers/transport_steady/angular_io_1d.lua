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
grid = meshgen:Execute(meshgen)

-- Set Material IDs
z_min = 0.0
z_max = widths[1]
for imat = 1, Nmat do
  z_max = z_min + widths[imat]
  log.Log(LOG_0, "imat=" .. imat .. ", zmin=" .. z_min .. ", zmax=" .. z_max)
  lv = logvol.RPPLogicalVolume.Create({ infx = true, infy = true, zmin = z_min, zmax = z_max })
  grid:SetBlockIDFromLogicalVolume(lv, imat - 1, true)
  z_min = z_max
end

-- Add cross sections to materials
total = { 50., 5., 0., 1., 1. }
c = { 0., 0., 0., 0.9, 0.9 }
xs_map = {}
for imat = 1, Nmat do
  xs_map[imat] = {
    block_ids = { imat - 1 },
    xs = xs.CreateSimpleOneGroup(total[imat], c[imat]),
  }
end

-- Create sources in 1st and 4th materials
src0 = lbs.VolumetricSource.Create({ block_ids = { 0 }, group_strength = { 50. } })
src1 = lbs.VolumetricSource.Create({ block_ids = { 3 }, group_strength = { 1. } })

-- Angular Quadrature
gl_quad = aquad.CreateGLProductQuadrature1DSlab(128)

-- LBS block option
num_groups = 1
lbs_block = {
  mesh = grid,
  num_groups = num_groups,
  groupsets = {
    {
      groups_from_to = { 0, num_groups - 1 },
      angular_quadrature = gl_quad,
      inner_linear_method = "petsc_gmres",
      l_abs_tol = 1.0e-9,
      l_max_its = 300,
      gmres_restart_interval = 30,
    },
  },
  xs_map = xs_map,
  options = {
    scattering_order = 0,
    spatial_discretization = "pwld",
    boundary_conditions = { { name = "zmin", type = "vacuum" }, { name = "zmax", type = "vacuum" } },
    save_angular_flux = true,
    volumetric_sources = { src0, src1 },
  },
}

phys1 = lbs.DiscreteOrdinatesProblem.Create(lbs_block)

-- Initialize and execute solver
ss_solver = lbs.SteadyStateSolver.Create({ lbs_problem = phys1 })

ss_solver:Initialize()
ss_solver:Execute()

leakage_left_1 = lbs.ComputeLeakage(phys1, { "zmin" })["zmin"][1]
leakage_right_1 = lbs.ComputeLeakage(phys1, { "zmax" })["zmax"][1]

lbs.WriteAngularFluxes(phys1, "angular_io")

phys2 = lbs.DiscreteOrdinatesProblem.Create(lbs_block)

ss_solver_2 = lbs.SteadyStateSolver.Create({ lbs_problem = phys2 })

ss_solver_2:Initialize()

lbs.ReadAngularFluxes(phys2, "angular_io")

leakage_left_2 = lbs.ComputeLeakage(phys2, { "zmin" })["zmin"][1]
leakage_right_2 = lbs.ComputeLeakage(phys2, { "zmax" })["zmax"][1]

leakage_left_diff = leakage_left_1 - leakage_left_2
leakage_right_diff = leakage_right_1 - leakage_right_2

log.Log(LOG_0, string.format("Leakage-Diff1=%.5e", leakage_left_diff))
log.Log(LOG_0, string.format("Leakage-Diff2=%.5e", leakage_right_diff))

os.execute("rm angular_io0.h5")
