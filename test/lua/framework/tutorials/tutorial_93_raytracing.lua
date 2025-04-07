-- Setup mesh
if nmesh == nil then
  nmesh = 11
end

nodes = {}
N = nmesh
L = 11.0
xmin = -L / 2
dx = L / N
for i = 1, (N + 1) do
  k = i - 1
  nodes[i] = xmin + k * dx
end

meshgen1 = mesh.OrthogonalMeshGenerator.Create({ node_sets = { nodes, nodes } })
grid = meshgen1:Execute()

-- Set block IDs
grid:SetUniformBlockID(0)

grid:SetOrthogonalBoundaries()

unit_sim_tests.SimTest93_RayTracing({ mesh = grid })

num_groups = 1
xs1g = xs.CreateSimpleOneGroup(1.0, 0.0)

---- Setup Physics
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
--pquad = aquad.CreateProductQuadrature(GAUSS_LEGENDRE_CHEBYSHEV,12*2*4, 12*4)
--aquad.OptimizeForPolarSymmetry(pquad, 4.0*math.pi)
--
----========== Groupset def
--gs0 = LBSCreateGroupset(phys1)
--cur_gs = gs0
--LBSGroupsetAddGroups(phys1,cur_gs,0,num_groups-1)
--LBSGroupsetSetQuadrature(phys1,cur_gs,pquad)
--LBSGroupsetSetAngleAggDiv(phys1,cur_gs,1)
--LBSGroupsetSetGroupSubsets(phys1,cur_gs,1)
--LBSGroupsetSetIterativeMethod(phys1,cur_gs,KRYLOV_RICHARDSON)
--LBSGroupsetSetResidualTolerance(phys1,cur_gs,1.0e-6)
--LBSGroupsetSetMaxIterations(phys1,cur_gs,0)
--LBSGroupsetSetGMRESRestartIntvl(phys1,cur_gs,100)
--
---- Set boundary conditions
--
---- Add point source
--src={}
--for g=1,num_groups do
--    src[g] = 0.0
--end
--src[1] = 1.0
--LBSAddPointSource(phys1, 0.0, 0.0, 0.0, src)
--
---- Set solver properties
--LBSSetProperty(phys1,DISCRETIZATION_METHOD,PWLD)
--LBSSetProperty(phys1,SCATTERING_ORDER,0)
--
---- Initialize and Execute Solver
--solver.Initialize(phys1)
--solver.Execute(phys1)

-- Add point source
src = {}
for g = 1, num_groups do
  src[g] = 0.0
end
src[1] = 1.0
pt_src = lbs.PointSource.Create({
  location = { 0.0, 0.0, 0.0 },
  strength = src,
})

-- Setup Physics
solver_name = "LBS"
pquad = aquad.CreateGLCProductQuadrature2DXY(12 * 4 * 2, 12 * 2 * 4 * 4)
lbs_block = {
  name = solver_name,
  mesh = grid,
  num_groups = num_groups,
  groupsets = {
    {
      groups_from_to = { 0, num_groups - 1 },
      angular_quadrature = pquad,
      inner_linear_method = "petsc_richardson",
      l_abs_tol = 1.0e-6,
      l_max_its = 0,
    },
  },
  xs_map = {
    { block_ids = { 0 }, xs = xs1g },
  },
  options = {
    scattering_order = 0,
    point_sources = { pt_src },
    field_function_prefix = solver_name,
  },
}

phys1 = lbs.DiscreteOrdinatesSolver.Create(lbs_block)

-- Initialize and Execute Solver
ss_solver = lbs.SteadyStateSolver.Create({ lbs_solver = phys1 })

ss_solver:Initialize()
ss_solver:Execute()

ff_m0 = lbs.GetScalarFieldFunctionList(phys1)

fieldfunc.ExportToVTKMulti({ ff_m0[1] }, "SimTest_93_LBS_" .. solver_name)
MPIBarrier()
if location_id == 0 then
  os.execute("rm SimTest_93*")
end
