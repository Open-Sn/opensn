-- 2D transport test in axialsymmetric cylindrical geometry with
-- vacuum boundary condition - monoenergetic.
-- SDM: PWLD
-- Test: Max-value=1.00000
num_procs = 4
--Structured mesh




--############################################### Check num_procs
if (check_num_procs==nil and number_of_processes ~= num_procs) then
    chiLog(LOG_0ERROR,"Incorrect amount of processors. " ..
                      "Expected "..tostring(num_procs)..
                      ". Pass check_num_procs=false to override if possible.")
    os.exit(false)
end

--############################################### Setup mesh
dim = 2
length = {1.0, 2.0, }
ncells = {50, 100, }
nodes = {}
for d = 1, dim do
  delta = length[d]/ncells[d]
  nodes[d] = {}
  for i = 0, ncells[d] do
    nodes[d][i+1] = i*delta
  end
end

meshgen1 = mesh.OrthogonalMeshGenerator.Create({ node_sets = {nodes[1],nodes[2]} })
mesh.MeshGenerator.Execute(meshgen1)

--############################################### Set Material IDs
vol0 = mesh.RPPLogicalVolume.Create
({ xmin=0.0,xmax=length[1],ymin=0.0,ymax=length[2], infz=true })
VolumeMesherSetProperty(MATID_FROMLOGICAL,vol0,0)

--############################################### Add materials
ngrp = 1
sigmat = 25.0
ratioc = 0.1
source = sigmat * (1 - ratioc)

material0 = PhysicsAddMaterial("Material_0");
PhysicsMaterialAddProperty(material0,TRANSPORT_XSECTIONS)
PhysicsMaterialAddProperty(material0,ISOTROPIC_MG_SOURCE)
PhysicsMaterialSetProperty(material0, TRANSPORT_XSECTIONS,
                           SIMPLEXS1, ngrp, sigmat, ratioc)
PhysicsMaterialSetProperty(material0, ISOTROPIC_MG_SOURCE,
                           SINGLE_VALUE, source)

--############################################### Setup Physics
pquad0 = chiCreateCylindricalProductQuadrature(GAUSS_LEGENDRE_CHEBYSHEV, 4, 8)

lbs_block =
{
  coord_system = 2,
  num_groups = ngrp,
  groupsets =
  {
    {
      groups_from_to = {0, ngrp-1},
      angular_quadrature_handle = pquad0,
      angle_aggregation_type = "azimuthal",
      inner_linear_method = "gmres",
      l_max_its = 100,
      l_abs_tol = 1.0e-12,
    }
  }
}

lbs_options =
{
  boundary_conditions = { { name = "xmin", type = "reflecting"} },
  scattering_order = 0,
}
phys1 = lbs.DiscreteOrdinatesCurvilinearSolver.Create(lbs_block)
lbs.SetOptions(phys1, lbs_options)

--############################################### Initialize and Execute Solver
ss_solver = lbs.SteadyStateSolver.Create({lbs_solver_handle = phys1})

SolverInitialize(ss_solver)
SolverExecute(ss_solver)


--phys0 = chiLBSCurvilinearCreateSolver(LBSCurvilinear.CYLINDRICAL)
--
----  angular quadrature
--pquad = chiCreateCylindricalProductQuadrature(GAUSS_LEGENDRE_CHEBYSHEV, 4, 8)
--
----  groups
--groups = {}
--for g = 1, ngrp do
--  groups[g] = LBSCreateGroup(phys0)
--end
--
----  groupsets
--gs0 = LBSCreateGroupset(phys0)
--LBSGroupsetAddGroups(phys0, gs0, 0, ngrp-1)
--LBSGroupsetSetQuadrature(phys0, gs0, pquad)
--LBSGroupsetSetAngleAggregationType(phys0, gs0, LBSGroupset.ANGLE_AGG_AZIMUTHAL)
--LBSGroupsetSetIterativeMethod(phys0, gs0, KRYLOV_GMRES_CYCLES)
--LBSGroupsetSetResidualTolerance(phys0, gs0, 1.0e-12)
--LBSGroupsetSetMaxIterations(phys0, gs0, 100)
--LBSGroupsetSetGMRESRestartIntvl(phys0, gs0, 30)
--
----  spatial discretisation
--LBSSetProperty(phys0, DISCRETIZATION_METHOD, PWLD)
--
----  scattering order
--LBSSetProperty(phys0, SCATTERING_ORDER, 0)
--
----------------------------------------------------------------------------------
----  boundary conditions
----------------------------------------------------------------------------------
--dirichlet_value = {}
--for g = 1, ngrp do
--  dirichlet_value[g] = 0
--end
--LBSSetProperty(phys0, BOUNDARY_CONDITION,
--                  XMIN, LBSBoundaryTypes.REFLECTING)
--LBSSetProperty(phys0, BOUNDARY_CONDITION,
--                  XMAX, LBSBoundaryTypes.INCIDENT_ISOTROPIC, dirichlet_value)
--LBSSetProperty(phys0, BOUNDARY_CONDITION,
--                  YMIN, LBSBoundaryTypes.INCIDENT_ISOTROPIC, dirichlet_value)
--LBSSetProperty(phys0, BOUNDARY_CONDITION,
--                  YMAX, LBSBoundaryTypes.INCIDENT_ISOTROPIC, dirichlet_value)
--
----------------------------------------------------------------------------------
----  solvers
----------------------------------------------------------------------------------
--SolverInitialize(phys0)
--SolverExecute(phys0)

--############################################### Exports
fflist, count = LBSGetScalarFieldFunctionList(phys1)
if master_export == nil then
  ExportMultiFieldFunctionToVTK(fflist, "ZRZPhi")
end

--############################################### Volume integrations
ffi1 = FFInterpolationCreate(VOLUME)
curffi = ffi1
FFInterpolationSetProperty(curffi, OPERATION, OP_MAX)
FFInterpolationSetProperty(curffi, LOGICAL_VOLUME, vol0)
FFInterpolationSetProperty(curffi, ADD_FIELDFUNCTION, fflist[1])

FFInterpolationInitialize(curffi)
FFInterpolationExecute(curffi)
maxval = FFInterpolationGetValue(curffi)

chiLog(LOG_0, string.format("Max-value=%.5f", maxval))
