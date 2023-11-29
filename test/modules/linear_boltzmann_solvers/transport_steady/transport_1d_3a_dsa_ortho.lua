-- 1D LinearBSolver test of a block of graphite with an air cavity. DSA and TG
-- SDM: PWLD
-- Test: WGS groups [0-62] Iteration    28 Residual 6.74299e-07 CONVERGED
-- and   WGS groups [63-167] Iteration    39 Residual 8.73816e-07 CONVERGED
num_procs = 4





--############################################### Check num_procs
if (check_num_procs==nil and number_of_processes ~= num_procs) then
  chiLog(LOG_0ERROR,"Incorrect amount of processors. " ..
    "Expected "..tostring(num_procs)..
    ". Pass check_num_procs=false to override if possible.")
  os.exit(false)
end

--############################################### Setup mesh
nodes={}
N=1000
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

meshgen1 = mesh.OrthogonalMeshGenerator.Create({ node_sets = {nodes} })
mesh.MeshGenerator.Execute(meshgen1)

--############################################### Set Material IDs
VolumeMesherSetMatIDToAll(0)

vol1 = mesh.RPPLogicalVolume.Create({infx=true, infy=true, zmin=-10.0, zmax=10.0})
VolumeMesherSetProperty(MATID_FROMLOGICAL,vol1,1)

--############################################### Add materials
materials = {}
materials[1] = PhysicsAddMaterial("Test Material");
materials[2] = PhysicsAddMaterial("Test Material2");

PhysicsMaterialAddProperty(materials[1],TRANSPORT_XSECTIONS)
PhysicsMaterialAddProperty(materials[2],TRANSPORT_XSECTIONS)

PhysicsMaterialAddProperty(materials[1],ISOTROPIC_MG_SOURCE)
PhysicsMaterialAddProperty(materials[2],ISOTROPIC_MG_SOURCE)


--num_groups = 1
--PhysicsMaterialSetProperty(materials[1],TRANSPORT_XSECTIONS,
--        SIMPLEXS1,num_groups,1.0,0.999)
num_groups = 168
PhysicsMaterialSetProperty(materials[1],TRANSPORT_XSECTIONS,
  CHI_XSFILE,"xs_graphite_pure.cxs")
PhysicsMaterialSetProperty(materials[2],TRANSPORT_XSECTIONS,
  CHI_XSFILE,"xs_air50RH.cxs")

src={}
for g=1,num_groups do
  src[g] = 0.0
end
src[1] = 1.0
PhysicsMaterialSetProperty(materials[1],ISOTROPIC_MG_SOURCE,FROM_ARRAY,src)
src[1] = 0.0
PhysicsMaterialSetProperty(materials[2],ISOTROPIC_MG_SOURCE,FROM_ARRAY,src)

--############################################### Setup Physics
pquad0 = CreateProductQuadrature(GAUSS_LEGENDRE_CHEBYSHEV,2, 2,false)
--pquad1 = CreateProductQuadrature(GAUSS_LEGENDRE_CHEBYSHEV,8, 8,false)

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
  scattering_order = 1,
}

phys1 = lbs.DiscreteOrdinatesSolver.Create(lbs_block)
lbs.SetOptions(phys1, lbs_options)

--############################################### Initialize and Execute Solver
ss_solver = lbs.SteadyStateSolver.Create({lbs_solver_handle = phys1})

SolverInitialize(ss_solver)
chiSolverExecute(ss_solver)

--############################################### Get field functions
fflist,count = chiLBSGetScalarFieldFunctionList(phys1)

--############################################### Exports
if (master_export == nil) then
  chiExportMultiFieldFunctionToVTK(fflist,"ZPhi")
end

--############################################### Plots
