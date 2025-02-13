-- 2D Transport test with localized material source
-- SDM: PWLD
-- Test: QoI Value=1.38399e-05
--       Inner Product=1.38405e-05
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

-- Set material IDs
grid:SetUniformMaterialID(0)

vol1a = logvol.RPPLogicalVolume.Create({
  infx = true,
  ymin = 0.0,
  ymax = 0.8 * L,
  infz = true,
})

grid:SetMaterialIDFromLogicalVolume(vol1a, 1, true)

vol0 = logvol.RPPLogicalVolume.Create({
  xmin = 2.5 - 0.166666,
  xmax = 2.5 + 0.166666,
  infy = true,
  infz = true,
})
grid:SetMaterialIDFromLogicalVolume(vol0, 0, true)

vol2 = logvol.RPPLogicalVolume.Create({
  xmin = 2.5 - 0.166666,
  xmax = 2.5 + 0.166666,
  ymin = 0.0,
  ymax = 2 * 0.166666,
  infz = true,
})
grid:SetMaterialIDFromLogicalVolume(vol2, 2, true)

vol1b = logvol.RPPLogicalVolume.Create({
  xmin = -1 + 2.5,
  xmax = 1 + 2.5,
  ymin = 0.9 * L,
  ymax = L,
  infz = true,
})
grid:SetMaterialIDFromLogicalVolume(vol1b, 1, true)

-- Create materials
materials = {}
materials[1] = mat.AddMaterial("Test Material1")
materials[2] = mat.AddMaterial("Test Material2")
materials[3] = mat.AddMaterial("Test Material3")

-- Add cross sections to materials
xs_1g1 = xs.CreateSimpleOneGroup(0.01, 0.01)
materials[1]:SetTransportXSections(xs_1g1)
xs_1g2 = xs.CreateSimpleOneGroup(0.1 * 20, 0.8)
materials[2]:SetTransportXSections(xs_1g2)
xs_1g3 = xs.CreateSimpleOneGroup(0.3 * 20, 0.0)
materials[3]:SetTransportXSections(xs_1g3)

-- Create sources
fwd_src = lbs.VolumetricSource.Create({ block_ids = { 2 }, group_strength = { 3.0 } })

-- Setup physics
pquad = aquad.CreateProductQuadrature(GAUSS_LEGENDRE_CHEBYSHEV, 48, 6)
aquad.OptimizeForPolarSymmetry(pquad, 4.0 * math.pi)

lbs_block = {
  mesh = grid,
  num_groups = 1,
  groupsets = {
    {
      groups_from_to = { 0, 0 },
      angular_quadrature = pquad,
      inner_linear_method = "petsc_gmres",
      l_abs_tol = 1.0e-6,
      l_max_its = 500,
      gmres_restart_interval = 100,
    },
  },
  options = {
    scattering_order = 0,
    volumetric_sources = { fwd_src },
  },
}
phys = lbs.DiscreteOrdinatesSolver.Create(lbs_block)

-- Forward solve
ss_solver = lbs.SteadyStateSolver.Create({ lbs_solver = phys })

ss_solver:Initialize()
ss_solver:Execute()

-- Get field functions
ff_m0 = fieldfunc.GetHandleByName("phi_g000_m00")
ff_m1 = fieldfunc.GetHandleByName("phi_g000_m01")
ff_m2 = fieldfunc.GetHandleByName("phi_g000_m02")

-- Define QoI region
qoi_vol = logvol.RPPLogicalVolume.Create({
  xmin = 0.5,
  xmax = 0.8333,
  ymin = 4.16666,
  ymax = 4.33333,
  infz = true,
})

-- Compute QoI
ffi = fieldfunc.FieldFunctionInterpolationVolume.Create()
ffi:SetOperationType(OP_SUM)
ffi:SetLogicalVolume(qoi_vol)
ffi:AddFieldFunction(ff_m0)

ffi:Initialize(ffi)
ffi:Execute(ffi)
fwd_qoi = ffi:GetValue()

-- Create adjoint source
adj_src = lbs.VolumetricSource.Create({ logical_volume = qoi_vol, group_strength = { 1.0 } })

-- Switch to adjoint mode
adjoint_options = {
  adjoint = true,
  volumetric_sources = { adj_src },
}
phys:SetOptions(adjoint_options)

-- Adjoint solve, write results
ss_solver:Execute()
lbs.WriteFluxMoments(phys, "adjoint_2d_1")

-- Create response evaluator
buffers = { { name = "buff", file_prefixes = { flux_moments = "adjoint_2d_1" } } }
mat_sources = { { material_id = 2, strength = { 3.0 } } }
response_options = {
  lbs_solver = phys,
  options = {
    buffers = buffers,
    sources = { material = mat_sources },
  },
}
evaluator = lbs.ResponseEvaluator.Create(response_options)

-- Evaluate response
adj_qoi = evaluator:EvaluateResponse("buff")

-- Print results
log.Log(LOG_0, string.format("QoI Value=%.5e", fwd_qoi))
log.Log(LOG_0, string.format("Inner Product=%.5e", adj_qoi))

-- Cleanup
MPIBarrier()
if location_id == 0 then
  os.execute("rm adjoint_2d_1*")
end
