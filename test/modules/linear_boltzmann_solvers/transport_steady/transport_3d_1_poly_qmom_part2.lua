-- 3D Transport test with Vacuum BCs and a material source reading source moments
-- SDM: PWLD
-- Test: Max-value=1.01701e-04 and 9.14681e-06
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
meshgen1 = mesh.ExtruderMeshGenerator.Create({
  inputs = {
    mesh.FromFileMeshGenerator.Create({
      filename = "../../../assets/mesh/SquareMesh2x2Quads.obj",
    }),
  },
  layers = { { z = 0.4, n = 2 }, { z = 0.8, n = 2 }, { z = 1.2, n = 2 }, { z = 1.6, n = 2 } }, -- layers
  partitioner = mesh.KBAGraphPartitioner.Create({
    nx = 2,
    ny = 2,
    xcuts = { 0.0 },
    ycuts = { 0.0 },
  }),
})
grid = meshgen1:Execute()

-- Set Material IDs
vol0 = logvol.RPPLogicalVolume.Create({ infx = true, infy = true, infz = true })
grid:SetMaterialIDFromLogicalVolume(vol0, 0, true)

vol1 = logvol.RPPLogicalVolume.Create({
  xmin = -0.5 / 8,
  xmax = 0.5 / 8,
  ymin = -0.5 / 8,
  ymax = 0.5 / 8,
  infz = true,
})
grid:SetMaterialIDFromLogicalVolume(vol1, 1, true)

-- Add materials
materials = {}
materials[1] = mat.AddMaterial("Test Material")
materials[2] = mat.AddMaterial("Test Material2")

num_groups = 21
xs_graphite = xs.LoadFromOpenSn("xs_graphite_pure.xs")
materials[1]:SetTransportXSections(xs_graphite)
materials[2]:SetTransportXSections(xs_graphite)

src = {}
for g = 1, num_groups do
  src[g] = 0.0
end

mg_src0 = xs.IsotropicMultiGroupSource.FromArray(src)
materials[1]:SetIsotropicMGSource(mg_src0)
src[1] = 1.0
mg_src1 = xs.IsotropicMultiGroupSource.FromArray(src)
materials[2]:SetIsotropicMGSource(mg_src1)

-- Setup Physics
pquad0 = aquad.CreateProductQuadrature(GAUSS_LEGENDRE_CHEBYSHEV, 2, 2)

lbs_block = {
  mesh = grid,
  num_groups = num_groups,
  groupsets = {
    {
      groups_from_to = { 0, 20 },
      angular_quadrature = pquad0,
      --angle_aggregation_type = "single",
      angle_aggregation_num_subsets = 1,
      groupset_num_subsets = 1,
      inner_linear_method = "petsc_gmres",
      l_abs_tol = 1.0e-6,
      l_max_its = 300,
      gmres_restart_interval = 100,
    },
  },
}

lbs_options = {
  scattering_order = 1,
  use_source_moments = true,
}

phys1 = lbs.DiscreteOrdinatesSolver.Create(lbs_block)
phys1:SetOptions(lbs_options)

-- Initialize and Execute Solver
ss_solver = lbs.SteadyStateSolver.Create({ lbs_solver = phys1 })
ss_solver:Initialize()
lbs.ReadSourceMoments(phys1, "Qmoms", false)

ss_solver:Execute()

-- Get field functions
fflist = lbs.GetScalarFieldFunctionList(phys1)

-- Slice plot
--slices = {}
--for k=1,count do
--    slices[k] = fieldfunc.FFInterpolationCreate(SLICE)
--    fieldfunc.SetProperty(slices[k],SLICE_POINT,{x = 0.0, y = 0.0, z = 0.8001})
--    fieldfunc.SetProperty(slices[k],ADD_FIELDFUNCTION,fflist[k])
--    --fieldfunc.SetProperty(slices[k],SLICE_TANGENT,{x = 0.393, y = 1.0-0.393, z = 0})
--    --fieldfunc.SetProperty(slices[k],SLICE_NORMAL,{x = -(1.0-0.393), y = -0.393, z = 0.0})
--    --fieldfunc.SetProperty(slices[k],SLICE_BINORM,{x = 0.0, y = 0.0, z = 1.0})
--    fieldfunc.Initialize(slices[k])
--    fieldfunc.Execute(slices[k])
--    fieldfunc.ExportToPython(slices[k])
--end

-- Volume integrations
ffi1 = fieldfunc.FieldFunctionInterpolationVolume.Create()
curffi = ffi1
curffi:SetOperationType(OP_MAX)
curffi:SetLogicalVolume(vol0)
curffi:AddFieldFunction(fflist[1])

curffi:Initialize()
curffi:Execute()
maxval = curffi:GetValue()

log.Log(LOG_0, string.format("Max-value1=%.5e", maxval))

ffi1 = fieldfunc.FieldFunctionInterpolationVolume.Create()
curffi = ffi1
curffi:SetOperationType(OP_MAX)
curffi:SetLogicalVolume(vol0)
curffi:AddFieldFunction(fflist[20])

curffi:Initialize()
curffi:Execute()
maxval = curffi:GetValue()

log.Log(LOG_0, string.format("Max-value2=%.5e", maxval))

-- Exports
if master_export == nil then
  fieldfunc.ExportToVTKMulti(fflist, "ZPhi3DColl")
end

-- Plots
if location_id == 0 and master_export == nil then
  --os.execute("python ZPFFI00.py")
  ----os.execute("python ZPFFI11.py")
  --local handle = io.popen("python ZPFFI00.py")
  print("Execution completed")
end
MPIBarrier()
if location_id == 0 then
  os.execute("rm Qmoms*")
end
