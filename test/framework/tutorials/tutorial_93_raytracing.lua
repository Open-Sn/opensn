--############################################### Setup mesh
if (nmesh==nil) then nmesh = 11 end

nodes={}
N=nmesh
L=11.0
xmin = -L/2
dx = L/N
for i=1,(N+1) do
    k=i-1
    nodes[i] = xmin + k*dx
end

meshgen1 = mesh.OrthogonalMeshGenerator.Create({ node_sets = {nodes,nodes} })
mesh.MeshGenerator.Execute(meshgen1)

--############################################### Set Material IDs
VolumeMesherSetMatIDToAll(0)

VolumeMesherSetupOrthogonalBoundaries()

unit_tests.SimTest93_RayTracing()


--###############################################
--############################################### Add materials
materials = {}
materials[1] = PhysicsAddMaterial("Test Material");

PhysicsMaterialAddProperty(materials[1],TRANSPORT_XSECTIONS)

num_groups = 1
PhysicsMaterialSetProperty(materials[1],
        TRANSPORT_XSECTIONS,
        SIMPLEXS0,1,0.27)



----############################################### Setup Physics
--solver_name = "LBS"
--phys1 = LBSCreateSolver(solver_name)
--
----========== Groups
--grp = {}
--for g=1,num_groups do
--    grp[g] = LBSCreateGroup(phys1)
--end
--
----========== ProdQuad
--pquad = CreateProductQuadrature(GAUSS_LEGENDRE_CHEBYSHEV,12*2*4, 12*4)
--OptimizeAngularQuadratureForPolarSymmetry(pquad, 4.0*math.pi)
--
----========== Groupset def
--gs0 = LBSCreateGroupset(phys1)
--cur_gs = gs0
--chiLBSGroupsetAddGroups(phys1,cur_gs,0,num_groups-1)
--chiLBSGroupsetSetQuadrature(phys1,cur_gs,pquad)
--chiLBSGroupsetSetAngleAggDiv(phys1,cur_gs,1)
--chiLBSGroupsetSetGroupSubsets(phys1,cur_gs,1)
--chiLBSGroupsetSetIterativeMethod(phys1,cur_gs,KRYLOV_RICHARDSON)
--chiLBSGroupsetSetResidualTolerance(phys1,cur_gs,1.0e-6)
--chiLBSGroupsetSetMaxIterations(phys1,cur_gs,0)
--chiLBSGroupsetSetGMRESRestartIntvl(phys1,cur_gs,100)
--
----############################################### Set boundary conditions
--
----############################################### Add point source
--src={}
--for g=1,num_groups do
--    src[g] = 0.0
--end
--src[1] = 1.0
--chiLBSAddPointSource(phys1, 0.0, 0.0, 0.0, src)
--
----############################################### Set solver properties
--chiLBSSetProperty(phys1,DISCRETIZATION_METHOD,PWLD)
--chiLBSSetProperty(phys1,SCATTERING_ORDER,0)
--
----############################################### Initialize and Execute Solver
--SolverInitialize(phys1)
--chiSolverExecute(phys1)






--############################################### Setup Physics
solver_name = "LBS"
pquad = CreateProductQuadrature(GAUSS_LEGENDRE_CHEBYSHEV,12*2*4, 12*4)
OptimizeAngularQuadratureForPolarSymmetry(pquad, 4.0*math.pi)
lbs_block =
{
    name = solver_name,
    num_groups = num_groups,
    groupsets =
    {
        {
            groups_from_to = {0, num_groups-1},
            angular_quadrature_handle = pquad,
            inner_linear_method = "richardson",
            l_abs_tol = 1.0e-6,
            l_max_its = 0,
        }
    },
    options = {scattering_order = 0, field_function_prefix = solver_name}
}

phys1 = lbs.DiscreteOrdinatesSolver.Create(lbs_block)

--############################################### Add point source
src={}
for g=1,num_groups do
    src[g] = 0.0
end
src[1] = 1.0
chiLBSAddPointSource(phys1, 0.0, 0.0, 0.0, src)

--############################################### Initialize and Execute Solver
ss_solver = lbs.SteadyStateSolver.Create({lbs_solver_handle = phys1})

SolverInitialize(ss_solver)
chiSolverExecute(ss_solver)

ff_m0 = chiLBSGetScalarFieldFunctionList(phys1)

chiExportMultiFieldFunctionToVTK({ff_m0[1]},"SimTest_93_LBS_"..solver_name)
MPIBarrier()
if (location_id == 0) then
    os.execute("rm SimTest_93*")
end
