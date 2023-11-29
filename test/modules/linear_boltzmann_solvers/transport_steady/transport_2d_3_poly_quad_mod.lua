-- 2D Transport test with Vacuum and Incident-isotropic BC. Quadrature optimized for polar symmetry.
-- SDM: PWLD
-- Test: Max-value=0.50758 and 2.52527e-04
num_procs = 4





--############################################### Check num_procs
if (check_num_procs==nil and number_of_processes ~= num_procs) then
  chiLog(LOG_0ERROR,"Incorrect amount of processors. " ..
    "Expected "..tostring(num_procs)..
    ". Pass check_num_procs=false to override if possible.")
  os.exit(false)
end

--############################################### Setup mesh
meshgen1 = mesh.MeshGenerator.Create
({
  inputs =
  {
    mesh.FromFileMeshGenerator.Create
    ({
      filename="../../../../resources/TestMeshes/SquareMesh2x2QuadsBlock.obj"
    }),
  },
  partitioner = chi.KBAGraphPartitioner.Create
  ({
    nx = 2, ny=2, nz=1,
    xcuts = {0.0}, ycuts = {0.0},
  })
})
mesh.MeshGenerator.Execute(meshgen1)

--############################################### Set Material IDs
chiVolumeMesherSetMatIDToAll(0)


--############################################### Add materials
materials = {}
materials[1] = chiPhysicsAddMaterial("Test Material");
materials[2] = chiPhysicsAddMaterial("Test Material2");

chiPhysicsMaterialAddProperty(materials[1],TRANSPORT_XSECTIONS)
chiPhysicsMaterialAddProperty(materials[2],TRANSPORT_XSECTIONS)

chiPhysicsMaterialAddProperty(materials[1],ISOTROPIC_MG_SOURCE)
chiPhysicsMaterialAddProperty(materials[2],ISOTROPIC_MG_SOURCE)


num_groups = 168
chiPhysicsMaterialSetProperty(materials[1],TRANSPORT_XSECTIONS,
  CHI_XSFILE,"xs_3_170.cxs")
chiPhysicsMaterialSetProperty(materials[2],TRANSPORT_XSECTIONS,
  CHI_XSFILE,"xs_3_170.cxs")

--chiPhysicsMaterialSetProperty(materials[1],TRANSPORT_XSECTIONS,SIMPLEXS0,num_groups,0.1)
--chiPhysicsMaterialSetProperty(materials[2],TRANSPORT_XSECTIONS,SIMPLEXS0,num_groups,0.1)

src={}
for g=1,num_groups do
  src[g] = 0.0
end
--src[1] = 1.0
chiPhysicsMaterialSetProperty(materials[1],ISOTROPIC_MG_SOURCE,FROM_ARRAY,src)
chiPhysicsMaterialSetProperty(materials[2],ISOTROPIC_MG_SOURCE,FROM_ARRAY,src)

--############################################### Setup Physics
pquad0 = chiCreateProductQuadrature(GAUSS_LEGENDRE_CHEBYSHEV,2, 1)
chiOptimizeAngularQuadratureForPolarSymmetry(pquad, 4.0*math.pi)

lbs_block =
{
  num_groups = num_groups,
  groupsets =
  {
    {
      groups_from_to = {0, 62},
      angular_quadrature_handle = pquad0,
      angle_aggregation_num_subsets = 1,
      groupset_num_subsets = 2,
      inner_linear_method = "gmres",
      l_abs_tol = 1.0e-6,
      l_max_its = 300,
      gmres_restart_interval = 100,
    },
    {
      groups_from_to = {63, num_groups-1},
      angular_quadrature_handle = pquad0,
      angle_aggregation_num_subsets = 1,
      groupset_num_subsets = 2,
      inner_linear_method = "gmres",
      l_abs_tol = 1.0e-6,
      l_max_its = 300,
      gmres_restart_interval = 100,
    },
  }
}
bsrc={}
for g=1,num_groups do
  bsrc[g] = 0.0
end
bsrc[1] = 1.0/4.0/math.pi

lbs_options =
{
  boundary_conditions =
  {
    {
      name = "xmin",
      type = "incident_isotropic",
      group_strength = bsrc
    }
  },
  scattering_order = 1,
}

phys1 = lbs.DiscreteOrdinatesSolver.Create(lbs_block)
lbs.SetOptions(phys1, lbs_options)

--############################################### Initialize and Execute Solver
ss_solver = lbs.SteadyStateSolver.Create({lbs_solver_handle = phys1})

chiSolverInitialize(ss_solver)
chiSolverExecute(ss_solver)

--############################################### Get field functions
fflist,count = chiLBSGetScalarFieldFunctionList(phys1)

--############################################### Slice plot
slice2 = chiFFInterpolationCreate(SLICE)
chiFFInterpolationSetProperty(slice2,SLICE_POINT,0.0,0.0,0.025)
chiFFInterpolationSetProperty(slice2,ADD_FIELDFUNCTION,fflist[1])

chiFFInterpolationInitialize(slice2)
chiFFInterpolationExecute(slice2)

--############################################### Volume integrations
vol0 = mesh.RPPLogicalVolume.Create({infx=true, infy=true, infz=true})

ffi1 = chiFFInterpolationCreate(VOLUME)
curffi = ffi1
chiFFInterpolationSetProperty(curffi,OPERATION,OP_MAX)
chiFFInterpolationSetProperty(curffi,LOGICAL_VOLUME,vol0)
chiFFInterpolationSetProperty(curffi,ADD_FIELDFUNCTION,fflist[1])

chiFFInterpolationInitialize(curffi)
chiFFInterpolationExecute(curffi)
maxval = chiFFInterpolationGetValue(curffi)

chiLog(LOG_0,string.format("Max-value1=%.5f", maxval))

--############################################### Volume integrations
ffi1 = chiFFInterpolationCreate(VOLUME)
curffi = ffi1
chiFFInterpolationSetProperty(curffi,OPERATION,OP_MAX)
chiFFInterpolationSetProperty(curffi,LOGICAL_VOLUME,vol0)
chiFFInterpolationSetProperty(curffi,ADD_FIELDFUNCTION,fflist[160])

chiFFInterpolationInitialize(curffi)
chiFFInterpolationExecute(curffi)
maxval = chiFFInterpolationGetValue(curffi)

chiLog(LOG_0,string.format("Max-value2=%.5e", maxval))

--############################################### Exports
if master_export == nil then
  chiFFInterpolationExportPython(slice2)
end

--############################################### Plots
if (location_id == 0 and master_export == nil) then
  local handle = io.popen("python ZPFFI00.py")
end
