-- 1D Transport test with Vacuum and Incident-isotropic BC.
-- SDM: PWLD
-- Test: Max-value=0.49903 and 7.18243e-4
num_procs = 3

--############################################### Check num_procs
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

--############################################### Setup mesh
nodes = {}
N = 100
L = 30.0
xmin = 0.0
dx = L / N
for i = 1, (N + 1) do
  k = i - 1
  nodes[i] = xmin + k * dx
end

meshgen1 = mesh.OrthogonalMeshGenerator.Create({ node_sets = { nodes } })
mesh.MeshGenerator.Execute(meshgen1)

--############################################### Set Material IDs
mesh.SetUniformMaterialID(0)

--############################################### Add materials
materials = {}
materials[1] = mat.AddMaterial("Test Material")
materials[2] = mat.AddMaterial("Test Material2")

num_groups = 168
mat.SetProperty(materials[1], TRANSPORT_XSECTIONS, OPENSN_XSFILE, "xs_3_170.xs")
mat.SetProperty(materials[2], TRANSPORT_XSECTIONS, OPENSN_XSFILE, "xs_3_170.xs")

src = {}
for g = 1, num_groups do
  src[g] = 0.0
end
--src[1] = 1.0
mat.SetProperty(materials[1], ISOTROPIC_MG_SOURCE, FROM_ARRAY, src)
mat.SetProperty(materials[2], ISOTROPIC_MG_SOURCE, FROM_ARRAY, src)

--############################################### Setup Physics
pquad0 = aquad.CreateProductQuadrature(GAUSS_LEGENDRE, 40)
lbs_block = {
  num_groups = num_groups,
  groupsets = {
    {
      groups_from_to = { 0, 62 },
      angular_quadrature_handle = pquad0,
      angle_aggregation_num_subsets = 1,
      groupset_num_subsets = 8,
      inner_linear_method = "petsc_gmres",
      l_abs_tol = 1.0e-6,
      l_max_its = 300,
      gmres_restart_interval = 100,
    },
    {
      groups_from_to = { 63, num_groups - 1 },
      angular_quadrature_handle = pquad0,
      angle_aggregation_num_subsets = 1,
      groupset_num_subsets = 8,
      inner_linear_method = "petsc_gmres",
      l_abs_tol = 1.0e-6,
      l_max_its = 300,
      gmres_restart_interval = 100,
    },
  },
}

bsrc = {}
for g = 1, num_groups do
  bsrc[g] = 0.0
end
bsrc[1] = 1.0 / 2

lbs_options = {
  boundary_conditions = {
    {
      name = "zmin",
      type = "isotropic",
      group_strength = bsrc,
    },
  },
  scattering_order = 5,
  max_ags_iterations = 1,
}

phys1 = lbs.DiscreteOrdinatesSolver.Create(lbs_block)
lbs.SetOptions(phys1, lbs_options)

--############################################### Initialize and Execute Solver
ss_solver = lbs.SteadyStateSolver.Create({ lbs_solver_handle = phys1 })

solver.Initialize(ss_solver)
solver.Execute(ss_solver)

--############################################### Get field functions
fflist, count = lbs.GetScalarFieldFunctionList(phys1)

--############################################### Line plot
--Testing consolidated interpolation
cline = fieldfunc.FFInterpolationCreate(LINE)
fieldfunc.SetProperty(cline, LINE_FIRSTPOINT, { x = 0.0, y = 0.0, z = 0.0001 + xmin })
fieldfunc.SetProperty(cline, LINE_SECONDPOINT, { x = 0.0, y = 0.0, z = 29.999 + xmin })
fieldfunc.SetProperty(cline, LINE_NUMBEROFPOINTS, 50)

for k = 165, 165 do
  fieldfunc.SetProperty(cline, ADD_FIELDFUNCTION, fflist[k])
end

fieldfunc.Initialize(cline)
fieldfunc.Execute(cline)

--############################################### Volume integrations
vol0 = logvol.RPPLogicalVolume.Create({ infx = true, infy = true, infz = true })
ffi1 = fieldfunc.FFInterpolationCreate(VOLUME)
curffi = ffi1
fieldfunc.SetProperty(curffi, OPERATION, OP_MAX)
fieldfunc.SetProperty(curffi, LOGICAL_VOLUME, vol0)
fieldfunc.SetProperty(curffi, ADD_FIELDFUNCTION, fflist[1])

fieldfunc.Initialize(curffi)
fieldfunc.Execute(curffi)
maxval = fieldfunc.GetValue(curffi)

log.Log(LOG_0, string.format("Max-value1=%.5f", maxval))

ffi2 = fieldfunc.FFInterpolationCreate(VOLUME)
curffi = ffi2
fieldfunc.SetProperty(curffi, OPERATION, OP_MAX)
fieldfunc.SetProperty(curffi, LOGICAL_VOLUME, vol0)
fieldfunc.SetProperty(curffi, ADD_FIELDFUNCTION, fflist[160])

fieldfunc.Initialize(curffi)
fieldfunc.Execute(curffi)
maxval = fieldfunc.GetValue(curffi)

log.Log(LOG_0, string.format("Max-value2=%.5e", maxval))

--############################################### Exports
if master_export == nil then
  fieldfunc.ExportToCSV(cline)
end

--############################################### Plots
if location_id == 0 and master_export == nil then
  local handle = io.popen("python3 ZLFFI00.py")
end
