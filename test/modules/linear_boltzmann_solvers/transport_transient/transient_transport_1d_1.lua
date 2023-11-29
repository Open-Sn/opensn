-- 1D Transient Transport test with Vacuum BC.
-- SDM: PWLD
-- Test:
num_procs = 2





--############################################### Check num_procs
if (check_num_procs==nil and number_of_processes ~= num_procs) then
    chiLog(LOG_0ERROR,"Incorrect amount of processors. " ..
                      "Expected "..tostring(num_procs)..
                      ". Pass check_num_procs=false to override if possible.")
    os.exit(false)
end

--############################################### Setup mesh
nodes={}
N=20
L=100.0
xmin = 0.0
dx = L/N
for i=1,(N+1) do
    k=i-1
    nodes[i] = xmin + k*dx
end

meshgen1 = mesh.OrthogonalMeshGenerator.Create({ node_sets = {nodes} })
mesh.MeshGenerator.Execute(meshgen1)

--############################################### Set Material IDs
VolumeMesherSetMatIDToAll(0)

--############################################### Add materials
materials = {}
materials[1] = PhysicsAddMaterial("Test Material");
materials[2] = PhysicsAddMaterial("Test Material2");

PhysicsMaterialAddProperty(materials[1],TRANSPORT_XSECTIONS)
PhysicsMaterialAddProperty(materials[2],TRANSPORT_XSECTIONS)

PhysicsMaterialAddProperty(materials[1],ISOTROPIC_MG_SOURCE)
PhysicsMaterialAddProperty(materials[2],ISOTROPIC_MG_SOURCE)

-- Define microscopic cross sections
xs = PhysicsTransportXSCreate()
xs_file = "tests/transport_transient/subcritical_1g.cxs"
PhysicsTransportXSSet(xs, CHI_XSFILE, xs_file)

-- Define macroscopic cross sections
--xs = PhysicsTransportXSMakeCombined({{xs, 0.00264086}}) -- just sub-critical
xs = PhysicsTransportXSMakeCombined({{xs, 0.0424459}}) -- just sub-critical

num_groups = 1
PhysicsMaterialSetProperty(materials[1],TRANSPORT_XSECTIONS,
        EXISTING,xs)
PhysicsMaterialSetProperty(materials[2],TRANSPORT_XSECTIONS,
        EXISTING,xs)

src={}
for g=1,num_groups do
    src[g] = 0.0
end
--src[1] = 1.0
PhysicsMaterialSetProperty(materials[1],ISOTROPIC_MG_SOURCE,FROM_ARRAY,src)
PhysicsMaterialSetProperty(materials[2],ISOTROPIC_MG_SOURCE,FROM_ARRAY,src)

--############################################### Setup Physics
phys1 = LBSCreateTransientSolver()

--========== Groups
grp = {}
for g=1,num_groups do
    grp[g] = LBSCreateGroup(phys1)
end

--========== ProdQuad
pquad = CreateProductQuadrature(GAUSS_LEGENDRE,16)

--========== Groupset def
gs0 = LBSCreateGroupset(phys1)
cur_gs = gs0
LBSGroupsetAddGroups(phys1,cur_gs,0,num_groups-1)
LBSGroupsetSetQuadrature(phys1,cur_gs,pquad)
LBSGroupsetSetAngleAggDiv(phys1,cur_gs,1)
LBSGroupsetSetGroupSubsets(phys1,cur_gs,8)
LBSGroupsetSetIterativeMethod(phys1,cur_gs,KRYLOV_GMRES)
LBSGroupsetSetResidualTolerance(phys1,cur_gs,1.0e-6)
LBSGroupsetSetMaxIterations(phys1,cur_gs,1000)
LBSGroupsetSetGMRESRestartIntvl(phys1,cur_gs,100)
--LBSGroupsetSetWGDSA(phys1,cur_gs,30,1.0e-4,false," ")
--LBSGroupsetSetTGDSA(phys1,cur_gs,30,1.0e-4,false," ")

--
----############################################### Set boundary conditions
--bsrc={}
--for g=1,num_groups do
--    bsrc[g] = 0.0
--end
--bsrc[1] = 1.0/2
----LBSSetProperty(phys1,BOUNDARY_CONDITION,ZMIN,LBSBoundaryTypes.INCIDENT_ISOTROPIC,bsrc);
--
LBSSetProperty(phys1,DISCRETIZATION_METHOD,PWLD)
LBSSetProperty(phys1,SCATTERING_ORDER,1)

LBKESSetProperty(phys1, "MAX_ITERATIONS", 100)
LBKESSetProperty(phys1, "TOLERANCE", 1.0e-8)

LBSSetProperty(phys1, USE_PRECURSORS, true)

--LBSSetProperty(phys1, VERBOSE_INNER_ITERATIONS, false)
LBSSetProperty(phys1, VERBOSE_INNER_ITERATIONS, false)
LBSSetProperty(phys1, VERBOSE_OUTER_ITERATIONS, true)


--############################################### Initialize and Execute Solver
SolverInitialize(phys1)

LBTSSetProperty(phys1, "TIMESTEP", 1e-1)
LBTSSetProperty(phys1, "VERBOSITY_LEVEL", 0)
LBTSSetProperty(phys1, "TIMESTEP_METHOD", "CRANK_NICHOLSON")

phys1name = SolverGetName(phys1);

--for k=1,2 do
--    --LBTSSetProperty(phys1, "INHIBIT_ADVANCE", true)
--    SolverStep(phys1)
--    FRf = chiLBSComputeFissionRate(phys1,"NEW")
--    FRi = chiLBSComputeFissionRate(phys1,"OLD")
--    dt = chiLBTSGetProperty(phys1, "TIMESTEP")
--    t = chiLBTSGetProperty(phys1, "TIME")
--    period = dt/math.log(FRf/FRi)
--    chiLog(LOG_0, string.format("%s time=%10.3g dt=%10.3g period=%10.3g", phys1name,t,dt,period))
--end

time = 0.0
time_stop = 20.0
k=0
while (time < time_stop) do
    k = k + 1
    SolverStep(phys1)
    FRf = chiLBSComputeFissionRate(phys1,"NEW")
    FRi = chiLBSComputeFissionRate(phys1,"OLD")
    dt = chiLBTSGetProperty(phys1, "TIMESTEP")
    time = chiLBTSGetProperty(phys1, "TIME")
    period = dt/math.log(FRf/FRi)
    chiLog(LOG_0, string.format("%s %4d time=%10.3g dt=%10.4g period=%10.3g FR=%10.3e",
            phys1name,k,time,dt,period,FRf))
end
