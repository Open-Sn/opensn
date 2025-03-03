-- 2D transport test in axialsymmetric cylindrical geometry with
-- vacuum boundary condition - multigroup with DSA.
-- SDM: PWLD
-- Test: Max-valueG1=1.00000, Max-valueG2=0.25000
num_procs = 4
--Structured mesh

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
dim = 2
length = { 1.0, 2.0 }
ncells = { 50, 100 }
nodes = {}
for d = 1, dim do
  delta = length[d] / ncells[d]
  nodes[d] = {}
  for i = 0, ncells[d] do
    nodes[d][i + 1] = i * delta
  end
end

meshgen1 = mesh.OrthogonalMeshGenerator.Create({ node_sets = { nodes[1], nodes[2] } })
grid = meshgen1:Execute()

-- Set Material IDs
vol0 = logvol.RPPLogicalVolume.Create({
  xmin = 0.0,
  xmax = length[1],
  ymin = 0.0,
  ymax = length[2],
  infz = true,
})
grid:SetMaterialIDFromLogicalVolume(vol0, 0, true)

-- Add materials
ngrp = 2
sigmat = 20.0
ratioc = 0.4
source = {}
source[1] = sigmat * (1 - 0.5 * ratioc)
for g = 2, ngrp do
  source[g] = 0.
end

material0 = mat.AddMaterial("Material_0")
xs_data = xs.LoadFromOpenSn("transport_2d_cyl_2_multigroup.xs")
material0:SetTransportXSections(xs_data)
mg_src = lbs.VolumetricSource.Create({ block_ids = { 0 }, group_strength = source })

-- Setup Physics
pquad0 = aquad.CreateGLCProductQuadrature2DRZ(4, 8)

lbs_block = {
  mesh = grid,
  coord_system = 2,
  num_groups = ngrp,
  groupsets = {
    {
      groups_from_to = { 0, ngrp - 1 },
      angular_quadrature = pquad0,
      angle_aggregation_type = "azimuthal",
      inner_linear_method = "petsc_gmres",
      l_max_its = 100,
      l_abs_tol = 1.0e-12,
      apply_wgdsa = true,
      wgdsa_l_abs_tol = 1.0e-9,
      wgdsa_l_max_its = 50,
    },
  },
}

lbs_options = {
  boundary_conditions = { { name = "xmin", type = "reflecting" } },
  scattering_order = 0,
  volumetric_sources = { mg_src },
}

phys1 = lbs.DiscreteOrdinatesCurvilinearSolver.Create(lbs_block)
phys1:SetOptions(lbs_options)

-- Initialize and Execute Solver
ss_solver = lbs.SteadyStateSolver.Create({ lbs_solver = phys1 })

ss_solver:Initialize()
ss_solver:Execute()

-- Exports
fflist = lbs.GetScalarFieldFunctionList(phys1)
if master_export == nil then
  fieldfunc.ExportToVTKMulti(fflist, "ZRZPhi")
end

--  volume integrations - energy group 1
ffi1 = fieldfunc.FieldFunctionInterpolationVolume.Create()
curffi = ffi1
curffi:SetOperationType(OP_MAX)
curffi:SetLogicalVolume(vol0)
curffi:AddFieldFunction(fflist[1])

curffi:Initialize()
curffi:Execute()
maxval = curffi:GetValue()

log.Log(LOG_0, string.format("Max-valueG1=%.5f", maxval))

--  volume integrations - energy group 2
ffi1 = fieldfunc.FieldFunctionInterpolationVolume.Create()
curffi = ffi1
curffi:SetOperationType(OP_MAX)
curffi:SetLogicalVolume(vol0)
curffi:AddFieldFunction(fflist[2])

curffi:Initialize()
curffi:Execute()
maxval = curffi:GetValue()

log.Log(LOG_0, string.format("Max-valueG2=%.5f", maxval))
