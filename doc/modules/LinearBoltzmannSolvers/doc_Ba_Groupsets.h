/** \defgroup LuaLBSGroupsets Groupsets

The code below is an example of a complete specification of a groupset.

\code
--===================================== Setup physics
phys1 = LBSCreateSolver()
chiSolverAddRegion(phys1,region1)

LBSSetProperty(phys1,DISCRETIZATION_METHOD,PWLD3D)
LBSSetProperty(phys1,SCATTERING_ORDER,1)

--========== Groups
grp = {}
for g=1,num_groups do
    grp[g] = LBSCreateGroup(phys1)
end

--========== ProdQuad
pquad0 = CreateProductQuadrature(GAUSS_LEGENDRE_CHEBYSHEV,2, 2)
pquad1 = CreateProductQuadrature(GAUSS_LEGENDRE_CHEBYSHEV,8, 8)

--========== Groupset def
gs0 = LBSCreateGroupset(phys1)

cur_gs = gs0
LBSGroupsetAddGroups(phys1,cur_gs,0,15)
LBSGroupsetSetQuadrature(phys1,cur_gs,pquad0)
LBSGroupsetSetAngleAggDiv(phys1,cur_gs,1)
LBSGroupsetSetGroupSubsets(phys1,cur_gs,1)
LBSGroupsetSetIterativeMethod(phys1,cur_gs,NPT_GMRES)
LBSGroupsetSetResidualTolerance(phys1,cur_gs,1.0e-4)
LBSGroupsetSetMaxIterations(phys1,cur_gs,300)
LBSGroupsetSetGMRESRestartIntvl(phys1,cur_gs,30)
LBSGroupsetSetWGDSA(phys1,cur_gs,30,1.0e-4,false," ")
LBSGroupsetSetTGDSA(phys1,cur_gs,30,1.0e-4,false," ")
\endcode

Groupsets segregate the code into pieces arranged by the number of groups
it contains. A great deal of care must be taken with intergroupset transfer
since the order in which the groupsets are executed determine what information
will be available to them.

\ingroup LBSUtilities*/
