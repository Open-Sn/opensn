-- 2D Transport test with Vacuum and Incident-isotropic BC.
-- SDM: PWLD
-- Test: Max-value1=3.18785
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

-- Setup mesh
nodes = {}
N = 40
L = 10.0
xmin = -L / 2
dx = L / N
for i = 1, (N + 1) do
  k = i - 1
  nodes[i] = xmin + k * dx
end

meshgen1 = mesh.OrthogonalMeshGenerator.Create({ node_sets = { nodes, nodes } })
grid = meshgen1:Execute()

-- Set Material IDs
grid:SetUniformMaterialID(0)
-- Add materials
materials = {}
materials[1] = mat.AddMaterial("Test Material")

num_groups = 1
xs_air = xs.LoadFromOpenSn("xs_air50RH.xs")
materials[1]:SetTransportXSections(xs_air)

strength = {}
for g = 1, num_groups do
  strength[g] = 0.0
end
--src[1] = 1.0
mg_src0 = lbs.VolumetricSource.Create({ block_ids = { 1 }, group_strength = strength })

-- Setup Physics
pquad0 = aquad.CreateProductQuadrature(GAUSS_LEGENDRE_CHEBYSHEV, 12, 2)
aquad.OptimizeForPolarSymmetry(pquad0, 4.0 * math.pi)

lbs_block = {
  mesh = grid,
  num_groups = num_groups,
  groupsets = {
    {
      groups_from_to = { 0, 0 },
      angular_quadrature = pquad0,
      angle_aggregation_num_subsets = 1,
      inner_linear_method = "petsc_gmres",
      l_abs_tol = 1.0e-6,
      l_max_its = 300,
      gmres_restart_interval = 100,
    },
  },
}

--int cell_global_id
--int material_id

--VecXYZ location (.x .y and .z)
--VecXYZ normal

--array<int>      quadrature_angle_indices
--array<VecXYZ>   quadrature_angle_vectors
--array<PhiTheta> quadrature_phi_theta_angles (PhiTheta.phi and PhiTheta.theta)
--array<int>      group_indices

--double          evaluation_time
function luaBoundaryFunctionA(
  cell_global_id,
  material_id,
  location,
  normal,
  quadrature_angle_indices,
  quadrature_angle_vectors,
  quadrature_phi_theta_angles,
  group_indices,
  time
)
  num_angles = rawlen(quadrature_angle_vectors)
  num_groups = rawlen(group_indices)
  psi = {}
  dof_count = 0

  for ni = 1, num_angles do
    omega = quadrature_angle_vectors[ni]
    phi_theta = quadrature_phi_theta_angles[ni]
    for gi = 1, num_groups do
      g = group_indices[gi]

      value = 1.0
      if location.y < 0.0 or omega.y < 0.0 then
        value = 0.0
      end

      dof_count = dof_count + 1
      psi[dof_count] = value
    end
  end

  return psi
end

lbs_options = {
  boundary_conditions = {
    {
      name = "xmin",
      type = "incident_anisotropic_heterogeneous",
      function_name = "luaBoundaryFunctionA",
    },
  },
  scattering_order = 1,
  volumetric_sources = { mg_src0 },
}

phys1 = lbs.DiscreteOrdinatesSolver.Create(lbs_block)
phys1:SetOptions(lbs_options)

-- Initialize and Execute Solver
ss_solver = lbs.SteadyStateSolver.Create({ lbs_solver = phys1 })

ss_solver:Initialize()
ss_solver:Execute()

-- Get field functions
fflist = lbs.GetScalarFieldFunctionList(phys1)

-- Slice plot
slice2 = fieldfunc.FFInterpolationCreate(SLICE)
fieldfunc.SetProperty(slice2, SLICE_POINT, { x = 0.0, y = 0.0, z = 0.025 })
fieldfunc.SetProperty(slice2, ADD_FIELDFUNCTION, fflist[1])

fieldfunc.Initialize(slice2)
fieldfunc.Execute(slice2)

---- Volume integrations
vol0 = logvol.RPPLogicalVolume.Create({ infx = true, infy = true, infz = true })
ffi1 = fieldfunc.FieldFunctionInterpolationVolume.Create()
curffi = ffi1
curffi:SetOperationType(OP_MAX)
curffi:SetLogicalVolume(vol0)
curffi:AddFieldFunction(fflist[1])

curffi:Initialize()
curffi:Execute()
maxval = curffi:GetValue()

log.Log(LOG_0, string.format("Max-value1=%.5f", maxval))

---- Volume integrations
--ffi1 = fieldfunc.FieldFunctionInterpolationVolume.Create()
--curffi = ffi1
--fieldfunc.SetProperty(curffi,OPERATION,OP_MAX)
--fieldfunc.SetProperty(curffi,LOGICAL_VOLUME,vol0)
--fieldfunc.SetProperty(curffi,ADD_FIELDFUNCTION,fflist[160])
--
--curffi:Initialize()
--curffi:Execute()
--maxval = curffi:GetValue()
--
--log.Log(LOG_0,string.format("Max-value2=%.5e", maxval))

-- Exports
if master_export == nil then
  fieldfunc.ExportToPython(slice2)
end

-- Plots
if location_id == 0 and master_export == nil then
  local handle = io.popen("python ZPFFI00.py")
end
