-- 2D Transport test with point source Multigroup FWD
-- SDM: PWLD
-- Test:
--  QoI Value[0]= 1.12687e-06
--  QoI Value[1]= 2.95934e-06
--  QoI Value[2]= 3.92975e-06
--  QoI Value[3]= 4.18474e-06
--  QoI Value[4]= 3.89649e-06
--  QoI Value[5]= 3.30482e-06
--  QoI Value[6]= 1.54506e-06
--  QoI Value[7]= 6.74868e-07
--  QoI Value[8]= 3.06178e-07
--  QoI Value[9]= 2.07284e-07
--  sum(QoI Value)= 2.21354e-05
--  Inner Product=3.30607e-06
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

-- Create mesh
N = 60
L = 5.0
ds = L / N

nodes = {}
for i = 0, N do
  nodes[i + 1] = i * ds
end
meshgen = mesh.OrthogonalMeshGenerator.Create({ node_sets = { nodes, nodes } })
grid = meshgen:Execute()

-- Set block IDs
grid:SetUniformBlockID(0)

vol1a = logvol.RPPLogicalVolume.Create({
  infx = true,
  ymin = 0.0,
  ymax = 0.8 * L,
  infz = true,
})

grid:SetBlockIDFromLogicalVolume(vol1a, 1, true)

vol0 = logvol.RPPLogicalVolume.Create({
  xmin = 2.5 - 0.166666,
  xmax = 2.5 + 0.166666,
  infy = true,
  infz = true,
})
grid:SetBlockIDFromLogicalVolume(vol0, 0, true)

vol1b = logvol.RPPLogicalVolume.Create({
  xmin = -1 + 2.5,
  xmax = 1 + 2.5,
  ymin = 0.9 * L,
  ymax = L,
  infz = true,
})
grid:SetBlockIDFromLogicalVolume(vol1b, 1, true)

-- Add cross sections to materials
num_groups = 10
xs_1 = xs.LoadFromOpenSn("response_2d_3_mat1.xs")
xs_2 = xs.LoadFromOpenSn("response_2d_3_mat2.xs")

-- Create sources
src = {}
for g = 1, num_groups do
  src[g] = 0.0
end
src[1] = 1.0

loc = { 1.25 - 0.5 * ds, 1.5 * ds, 0.0 }
pt_src = lbs.PointSource.Create({ location = loc, strength = src })

-- Setup physics
pquad = aquad.CreateGLCProductQuadrature2DXY(4, 48)

lbs_block = {
  mesh = grid,
  num_groups = num_groups,
  groupsets = {
    {
      groups_from_to = { 0, num_groups - 1 },
      angular_quadrature = pquad,
      inner_linear_method = "petsc_gmres",
      l_abs_tol = 1.0e-6,
      l_max_its = 500,
      gmres_restart_interval = 100,
    },
  },
  xs_map = {
    { block_ids = { 0 }, xs = xs_1 },
    { block_ids = { 1 }, xs = xs_2 },
  },
  options = {
    scattering_order = 0,
    point_sources = { pt_src },
  },
}
phys = lbs.DiscreteOrdinatesSolver.Create(lbs_block)

-- Forward solve
ss_solver = lbs.SteadyStateSolver.Create({ lbs_solver = phys })

ss_solver:Initialize()
ss_solver:Execute()

-- Define QoI region
qoi_vol = logvol.RPPLogicalVolume.Create({
  xmin = 0.5,
  xmax = 0.8333,
  ymin = 4.16666,
  ymax = 4.33333,
  infz = true,
})

-- Compute QoI
fwd_qois = {}
fwd_qoi_sum = 0.0
for g = 0, num_groups - 1 do
  ff = fieldfunc.GetHandleByName(
    "phi_g" .. string.format("%03d", g) .. "_m" .. string.format("%02d", 0)
  )
  ffi = fieldfunc.FieldFunctionInterpolationVolume()
  ffi:SetOperationType(OP_SUM)
  ffi:SetLogicalVolume(qoi_vol)
  ffi:AddFieldFunction(ff)

  ffi:Initialize()
  ffi:Execute()
  fwd_qois[g + 1] = ffi:GetValue()

  fwd_qoi_sum = fwd_qoi_sum + fwd_qois[g + 1]
end

-- Create adjoint source
function ResponseFunction(xyz, mat_id)
  response = {}
  for g = 1, num_groups do
    if g == 6 then
      response[g] = 1.0
    else
      response[g] = 0.0
    end
  end
  return response
end

response_func = LuaVectorSpatialFunction.Create({ function_name = "ResponseFunction" })

adjoint_source = lbs.VolumetricSource.Create({
  logical_volume = qoi_vol,
  func = response_func,
})

-- Switch to adjoint mode
adjoint_options = {
  adjoint = true,
  volumetric_sources = { adjoint_source },
}
phys:SetOptions(adjoint_options)

-- Adjoint solve, write results
ss_solver:Execute()
lbs.WriteFluxMoments(phys, "adjoint_2d_3")

-- Create response evaluator
buffers = { { name = "buff", file_prefixes = { flux_moments = "adjoint_2d_3" } } }
pt_sources = { pt_src }
response_options = {
  lbs_solver = phys,
  options = {
    buffers = buffers,
    sources = { point = pt_sources },
  },
}
evaluator = lbs.ResponseEvaluator.Create(response_options)

-- Evaluate response
response = evaluator:EvaluateResponse("buff")

-- Print results
for g = 1, num_groups do
  pref = "QoI Value[" .. tostring(g - 1) .. "]"
  log.Log(LOG_0, string.format(pref .. "= %.5e", fwd_qois[g]))
end
log.Log(LOG_0, string.format("sum(QoI Values)= %.5e", fwd_qoi_sum))
log.Log(LOG_0, string.format("Inner Product=%.5e", response))

-- Cleanup
MPIBarrier()
if location_id == 0 then
  os.execute("rm adjoint_2d_3*")
end
