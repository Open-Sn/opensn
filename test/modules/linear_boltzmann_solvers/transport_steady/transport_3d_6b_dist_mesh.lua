-- 3D Transport test with distributedt-mesh + 2D ortho mesh + extruded mesh.
-- SDM: PWLD
-- Test: Max-value1=6.55387e+00
--       Max-value2=1.02940e+00

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

-- Cells
div = 8
Nx = math.floor(128 / div)
Ny = math.floor(128 / div)
Nz = math.floor(256 / div)

-- Dimensions
Lx = 10.0
Ly = 10.0
Lz = 10.0

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

meshgen1 = mesh.DistributedMeshGenerator.Create({
  inputs = {
    mesh.OrthogonalMeshGenerator.Create({ node_sets = { xmesh, ymesh } }),
    mesh.ExtruderMeshGenerator.Create({
      layers = { { z = Lz, n = Nz } },
    }),
  },
})

grid = meshgen1:Execute()

grid:SetUniformMaterialID(0)

-- Add materials
materials = {}
materials[1] = mat.AddMaterial("Test Material")

num_groups = 21
xs_graphite = xs.LoadFromOpenSn("xs_graphite_pure.xs")
materials[1]:SetTransportXSections(xs_graphite)

strength = {}
for g = 1, num_groups do
  strength[g] = 0.0
end
mg_src = lbs.VolumetricSource.Create({ block_ids = { 1 }, group_strength = strength })

-- Setup Physics
pquad0 = aquad.CreateProductQuadrature(GAUSS_LEGENDRE_CHEBYSHEV, 2, 4, false)

lbs_block = {
  mesh = grid,
  num_groups = num_groups,
  groupsets = {
    {
      groups_from_to = { 0, 20 },
      angular_quadrature = pquad0,
      angle_aggregation_type = "polar",
      angle_aggregation_num_subsets = 1,
      groupset_num_subsets = 1,
      inner_linear_method = "petsc_gmres",
      l_abs_tol = 1.0e-6,
      l_max_its = 300,
      gmres_restart_interval = 100,
    },
  },
  sweep_type = "CBC",
}
bsrc = {}
for g = 1, num_groups do
  bsrc[g] = 0.0
end
bsrc[1] = 1.0 / 4.0 / math.pi
lbs_options = {
  boundary_conditions = {
    { name = "xmin", type = "isotropic", group_strength = bsrc },
  },
  scattering_order = 1,
  save_angular_flux = true,
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

pp1 = post.CellVolumeIntegralPostProcessor.Create({
  name = "max-grp0",
  field_function = fflist[1],
  compute_volume_average = true,
  print_numeric_format = "scientific",
})
pp2 = post.CellVolumeIntegralPostProcessor.Create({
  name = "max-grp19",
  field_function = fflist[20],
  compute_volume_average = true,
  print_numeric_format = "scientific",
})
post.Execute({ pp1, pp2 })

if master_export == nil then
  fieldfunc.ExportToVTKMulti(fflist, "ZPhi")
end
