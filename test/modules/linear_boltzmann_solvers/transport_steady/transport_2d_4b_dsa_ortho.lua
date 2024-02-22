-- 2D LinearBSolver Same as 4a but with reflective BCs. DSA and TG
-- SDM: PWLD
-- Test: WGS groups [0-62] Iteration    54 Residual 5.00021e-07 CONVERGED
-- and   WGS groups [63-167] Iteration    56 Residual 9.73954e-07 CONVERGED
num_procs = 4





--############################################### Check num_procs
if (check_num_procs==nil and number_of_processes ~= num_procs) then
  log.Log(LOG_0ERROR,"Incorrect amount of processors. " ..
    "Expected "..tostring(num_procs)..
    ". Pass check_num_procs=false to override if possible.")
  os.exit(false)
end

--############################################### Setup mesh
nodes={}
N=20
L=100
--N=10
--L=200e6
xmin = -L/2
--xmin = 0.0
dx = L/N
for i=1,(N+1) do
  k=i-1
  nodes[i] = xmin + k*dx
end

meshgen1 = mesh.OrthogonalMeshGenerator.Create({ node_sets = {nodes,nodes} })
mesh.MeshGenerator.Execute(meshgen1)

--############################################### Set Material IDs
mesh.SetUniformMaterialID(0)

vol1 = mesh.RPPLogicalVolume.Create
({ xmin=-10.0,xmax=10.0,ymin=-10.0,ymax=10.0, infz=true })
mesh.SetMaterialIDFromLogicalVolume(vol1,1)

--############################################### Add materials
materials = {}
materials[1] = mat.AddMaterial("Test Material");
materials[2] = mat.AddMaterial("Test Material2");

mat.AddProperty(materials[1], TRANSPORT_XSECTIONS)
mat.AddProperty(materials[2], TRANSPORT_XSECTIONS)

mat.AddProperty(materials[1], ISOTROPIC_MG_SOURCE)
mat.AddProperty(materials[2], ISOTROPIC_MG_SOURCE)


--num_groups = 1
--mat.SetProperty(materials[1], TRANSPORT_XSECTIONS, SIMPLEXS1, num_groups, 1.0, 0.999)
num_groups = 168
mat.SetProperty(materials[1], TRANSPORT_XSECTIONS, OPENSN_XSFILE, "xs_graphite_pure.xs")
mat.SetProperty(materials[2], TRANSPORT_XSECTIONS, OPENSN_XSFILE, "xs_air50RH.xs")

src={}
for g=1,num_groups do
  src[g] = 0.0
end
src[1] = 1.0
mat.SetProperty(materials[1], ISOTROPIC_MG_SOURCE, FROM_ARRAY, src)
src[1] = 0.0
mat.SetProperty(materials[2], ISOTROPIC_MG_SOURCE, FROM_ARRAY, src)

--############################################### Setup Physics
pquad0 = aquad.CreateProductQuadrature(GAUSS_LEGENDRE_CHEBYSHEV, 2, 2, false)
aquad.OptimizeForPolarSymmetry(pquad0, 4.0*math.pi)

lbs_block =
{
  num_groups = num_groups,
  groupsets =
  {
    {
      groups_from_to = {0, 62},
      angular_quadrature_handle = pquad0,
      angle_aggregation_num_subsets = 1,
      groupset_num_subsets = 1,
      inner_linear_method = "gmres",
      l_abs_tol = 1.0e-6,
      l_max_its = 1000,
      gmres_restart_interval = 30,
      apply_wgdsa = true,
      wgdsa_l_abs_tol = 1.0e-2,
    },
    {
      groups_from_to = {63, num_groups-1},
      angular_quadrature_handle = pquad0,
      angle_aggregation_num_subsets = 1,
      groupset_num_subsets = 1,
      inner_linear_method = "gmres",
      l_abs_tol = 1.0e-6,
      l_max_its = 1000,
      gmres_restart_interval = 30,
      apply_wgdsa = true,
      apply_tgdsa = true,
      wgdsa_l_abs_tol = 1.0e-2,
    },
  }
}

lbs_options =
{
  boundary_conditions =
  {
   {name = "xmin",type = "reflecting"},
   {name = "ymin",type = "reflecting"},
  },
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

--############################################### Exports
if (master_export == nil) then
  fieldfunc.ExportToVTKMulti(fflist,"ZPhi")
end

--############################################### Plots
