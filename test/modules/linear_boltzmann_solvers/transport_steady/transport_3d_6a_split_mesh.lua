-- 3D Transport test with split-mesh + ortho mesh.
-- SDM: PWLD
-- Test: max-grp0(latest) =  1.131566e-01
--       max-grp19(latest) = 7.340585e-04

num_procs = 4





--############################################### Check num_procs
if (check_num_procs==nil and number_of_processes ~= num_procs) then
  Log(LOG_0ERROR,"Incorrect amount of processors. " ..
    "Expected "..tostring(num_procs)..
    ". Pass check_num_procs=false to override if possible.")
  os.exit(false)
end

-- Cells
div = 8
Nx = math.floor(128/div)
Ny = math.floor(128/div)
Nz = math.floor(256/div)

-- Dimensions
Lx = 10.0
Ly = 10.0
Lz = 10.0

xmesh = {}
xmin = 0.0
dx = Lx/Nx
for i = 1, (Nx+1) do
  k = i-1
  xmesh[i] = xmin + k*dx
end

ymesh = {}
ymin = 0.0
dy = Ly/Ny
for i = 1, (Ny+1) do
  k = i-1
  ymesh[i] = ymin + k*dy
end

zmesh = {}
zmin = 0.0
dz = Lz/Nz
for i = 1, (Nz+1) do
  k = i-1
  zmesh[i] = zmin + k*dz
end

meshgen1 = mesh.SplitFileMeshGenerator.Create
({
  inputs =
  {
    mesh.OrthogonalMeshGenerator.Create({ node_sets = {xmesh,ymesh,zmesh} })
  },
})

mesh.MeshGenerator.Execute(meshgen1)

mesh.SetMatIDToAll(0)

--############################################### Add materials
materials = {}
materials[1] = PhysicsAddMaterial("Test Material");

PhysicsMaterialAddProperty(materials[1],TRANSPORT_XSECTIONS)

PhysicsMaterialAddProperty(materials[1],ISOTROPIC_MG_SOURCE)


num_groups = 21
PhysicsMaterialSetProperty(materials[1],TRANSPORT_XSECTIONS,
  OPENSN_XSFILE,"xs_graphite_pure.xs")

src={}
for g=1,num_groups do
  src[g] = 0.0
end
PhysicsMaterialSetProperty(materials[1],ISOTROPIC_MG_SOURCE,FROM_ARRAY,src)

--############################################### Setup Physics
pquad0 = CreateProductQuadrature(GAUSS_LEGENDRE_CHEBYSHEV,2, 4)

lbs_block =
{
  num_groups = num_groups,
  groupsets =
  {
    {
      groups_from_to = {0, 20},
      angular_quadrature_handle = pquad0,
      angle_aggregation_type = "polar",
      angle_aggregation_num_subsets = 1,
      groupset_num_subsets = 1,
      inner_linear_method = "gmres",
      l_abs_tol = 1.0e-6,
      l_max_its = 300,
      gmres_restart_interval = 100,
    },
  },
  sweep_type = "CBC",
}
bsrc={}
for g=1,num_groups do
  bsrc[g] = 0.0
end
bsrc[1] = 1.0/4.0/math.pi;
lbs_options =
{
  boundary_conditions = { { name = "xmin", type = "incident_isotropic",
                            group_strength=bsrc}},
  scattering_order = 1,
  save_angular_flux = true
}

phys1 = lbs.DiscreteOrdinatesSolver.Create(lbs_block)
lbs.SetOptions(phys1, lbs_options)

--############################################### Initialize and Execute Solver
ss_solver = lbs.SteadyStateSolver.Create({lbs_solver_handle = phys1})

SolverInitialize(ss_solver)
SolverExecute(ss_solver)

--############################################### Get field functions
fflist,count = LBSGetScalarFieldFunctionList(phys1)

pp1 = CellVolumeIntegralPostProcessor.Create
({
  name="max-grp0",
  field_function = fflist[1],
  compute_volume_average = true,
  print_numeric_format = "scientific"
})
pp2 = CellVolumeIntegralPostProcessor.Create
({
  name="max-grp19",
  field_function = fflist[20],
  compute_volume_average = true,
  print_numeric_format = "scientific"
})
ExecutePostProcessors({ pp1, pp2 })

if (master_export == nil) then
  ExportMultiFieldFunctionToVTK(fflist,"ZPhi")
end

LogPrintTimingGraph()
