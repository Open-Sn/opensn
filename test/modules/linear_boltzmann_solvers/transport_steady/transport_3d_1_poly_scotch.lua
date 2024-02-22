-- 3D Transport test with Vacuum and Incident-isotropic BC.
-- SDM: PWLD
-- Test: Max-value=0.49903 and 7.18243e-4
num_procs = 4





--############################################### Check num_procs
if (check_num_procs==nil and number_of_processes ~= num_procs) then
  log.Log(LOG_0ERROR,"Incorrect amount of processors. " ..
    "Expected "..tostring(num_procs)..
    ". Pass check_num_procs=false to override if possible.")
  os.exit(false)
end

--############################################### Setup mesh
Nxy = 32
nodesxy = {}
dxy = 2/Nxy
dz = 1.6/8
for i=0,(Nxy) do
  nodesxy[i+1] = -1.0 + i*dxy
end
nodesz = {}
for k=0,8 do
  nodesz[k+1] = 0.0 + k*dz
end

meshgen = mesh.MeshGenerator.Create
({
  inputs =
  {
    mesh.OrthogonalMeshGenerator.Create
    ({
      node_sets = {nodesxy, nodesxy, nodesz}
    }),
  },
  partitioner = PETScGraphPartitioner.Create({type="ptscotch"})
})
mesh.MeshGenerator.Execute(meshgen)

--############################################### Set Material IDs
vol0 = mesh.RPPLogicalVolume.Create({infx=true, infy=true, infz=true})
mesh.SetMaterialIDFromLogicalVolume(vol0,0)

vol1 = mesh.RPPLogicalVolume.Create
({ xmin=-0.5,xmax=0.5,ymin=-0.5,ymax=0.5, infz=true })
mesh.SetMaterialIDFromLogicalVolume(vol1,1)


--############################################### Add materials
materials = {}
materials[1] = mat.AddMaterial("Test Material");
materials[2] = mat.AddMaterial("Test Material2");

mat.AddProperty(materials[1],TRANSPORT_XSECTIONS)
mat.AddProperty(materials[2],TRANSPORT_XSECTIONS)

mat.AddProperty(materials[1],ISOTROPIC_MG_SOURCE)
mat.AddProperty(materials[2],ISOTROPIC_MG_SOURCE)


num_groups = 21
mat.SetProperty(materials[1],TRANSPORT_XSECTIONS,
  OPENSN_XSFILE,"xs_graphite_pure.xs")
mat.SetProperty(materials[2],TRANSPORT_XSECTIONS,
  OPENSN_XSFILE,"xs_graphite_pure.xs")

src={}
for g=1,num_groups do
  src[g] = 0.0
end

mat.SetProperty(materials[1],ISOTROPIC_MG_SOURCE,FROM_ARRAY,src)
mat.SetProperty(materials[2],ISOTROPIC_MG_SOURCE,FROM_ARRAY,src)



--############################################### Setup Physics
pquad0 = aquad.CreateProductQuadrature(GAUSS_LEGENDRE_CHEBYSHEV,2, 2)

lbs_block =
{
  num_groups = num_groups,
  groupsets =
  {
    {
      groups_from_to = {0, 20},
      angular_quadrature_handle = pquad0,
      angle_aggregation_type = "single",
      angle_aggregation_num_subsets = 1,
      groupset_num_subsets = 1,
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
bsrc[1] = 1.0/4.0/math.pi;
lbs_options =
{
  boundary_conditions = { { name = "zmin", type = "isotropic",
                            group_strength=bsrc}},
  scattering_order = 1,
}

phys1 = lbs.DiscreteOrdinatesSolver.Create(lbs_block)
lbs.SetOptions(phys1, lbs_options)

--############################################### Initialize and Execute Solver
ss_solver = lbs.SteadyStateSolver.Create({lbs_solver_handle = phys1})

solver.Initialize(ss_solver)
solver.Execute(ss_solver)

--############################################### Get field functions
fflist,count = lbs.GetScalarFieldFunctionList(phys1)

--############################################### Volume integrations
ffi1 = fieldfunc.FFInterpolationCreate(VOLUME)
curffi = ffi1
fieldfunc.SetProperty(curffi,OPERATION,OP_MAX)
fieldfunc.SetProperty(curffi,LOGICAL_VOLUME,vol0)
fieldfunc.SetProperty(curffi,ADD_FIELDFUNCTION,fflist[1])

fieldfunc.Initialize(curffi)
fieldfunc.Execute(curffi)
maxval = fieldfunc.GetValue(curffi)

log.Log(LOG_0,string.format("Max-value1=%.5e", maxval))

ffi1 = fieldfunc.FFInterpolationCreate(VOLUME)
curffi = ffi1
fieldfunc.SetProperty(curffi,OPERATION,OP_MAX)
fieldfunc.SetProperty(curffi,LOGICAL_VOLUME,vol0)
fieldfunc.SetProperty(curffi,ADD_FIELDFUNCTION,fflist[20])

fieldfunc.Initialize(curffi)
fieldfunc.Execute(curffi)
maxval = fieldfunc.GetValue(curffi)

log.Log(LOG_0,string.format("Max-value2=%.5e", maxval))
