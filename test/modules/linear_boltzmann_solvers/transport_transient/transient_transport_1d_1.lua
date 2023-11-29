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
chiVolumeMesherSetMatIDToAll(0)

--############################################### Add materials
materials = {}
materials[1] = chiPhysicsAddMaterial("Test Material");
materials[2] = chiPhysicsAddMaterial("Test Material2");

chiPhysicsMaterialAddProperty(materials[1],TRANSPORT_XSECTIONS)
chiPhysicsMaterialAddProperty(materials[2],TRANSPORT_XSECTIONS)

chiPhysicsMaterialAddProperty(materials[1],ISOTROPIC_MG_SOURCE)
chiPhysicsMaterialAddProperty(materials[2],ISOTROPIC_MG_SOURCE)

-- Define microscopic cross sections
xs = chiPhysicsTransportXSCreate()
xs_file = "tests/transport_transient/subcritical_1g.cxs"
chiPhysicsTransportXSSet(xs, CHI_XSFILE, xs_file)

-- Define macroscopic cross sections
--xs = chiPhysicsTransportXSMakeCombined({{xs, 0.00264086}}) -- just sub-critical
xs = chiPhysicsTransportXSMakeCombined({{xs, 0.0424459}}) -- just sub-critical

num_groups = 1
chiPhysicsMaterialSetProperty(materials[1],TRANSPORT_XSECTIONS,
        EXISTING,xs)
chiPhysicsMaterialSetProperty(materials[2],TRANSPORT_XSECTIONS,
        EXISTING,xs)

src={}
for g=1,num_groups do
    src[g] = 0.0
end
--src[1] = 1.0
chiPhysicsMaterialSetProperty(materials[1],ISOTROPIC_MG_SOURCE,FROM_ARRAY,src)
chiPhysicsMaterialSetProperty(materials[2],ISOTROPIC_MG_SOURCE,FROM_ARRAY,src)

--############################################### Setup Physics
phys1 = chiLBSCreateTransientSolver()

--========== Groups
grp = {}
for g=1,num_groups do
    grp[g] = chiLBSCreateGroup(phys1)
end

--========== ProdQuad
pquad = chiCreateProductQuadrature(GAUSS_LEGENDRE,16)

--========== Groupset def
gs0 = chiLBSCreateGroupset(phys1)
cur_gs = gs0
chiLBSGroupsetAddGroups(phys1,cur_gs,0,num_groups-1)
chiLBSGroupsetSetQuadrature(phys1,cur_gs,pquad)
chiLBSGroupsetSetAngleAggDiv(phys1,cur_gs,1)
chiLBSGroupsetSetGroupSubsets(phys1,cur_gs,8)
chiLBSGroupsetSetIterativeMethod(phys1,cur_gs,KRYLOV_GMRES)
chiLBSGroupsetSetResidualTolerance(phys1,cur_gs,1.0e-6)
chiLBSGroupsetSetMaxIterations(phys1,cur_gs,1000)
chiLBSGroupsetSetGMRESRestartIntvl(phys1,cur_gs,100)
--chiLBSGroupsetSetWGDSA(phys1,cur_gs,30,1.0e-4,false," ")
--chiLBSGroupsetSetTGDSA(phys1,cur_gs,30,1.0e-4,false," ")

--
----############################################### Set boundary conditions
--bsrc={}
--for g=1,num_groups do
--    bsrc[g] = 0.0
--end
--bsrc[1] = 1.0/2
----chiLBSSetProperty(phys1,BOUNDARY_CONDITION,ZMIN,LBSBoundaryTypes.INCIDENT_ISOTROPIC,bsrc);
--
chiLBSSetProperty(phys1,DISCRETIZATION_METHOD,PWLD)
chiLBSSetProperty(phys1,SCATTERING_ORDER,1)

chiLBKESSetProperty(phys1, "MAX_ITERATIONS", 100)
chiLBKESSetProperty(phys1, "TOLERANCE", 1.0e-8)

chiLBSSetProperty(phys1, USE_PRECURSORS, true)

--chiLBSSetProperty(phys1, VERBOSE_INNER_ITERATIONS, false)
chiLBSSetProperty(phys1, VERBOSE_INNER_ITERATIONS, false)
chiLBSSetProperty(phys1, VERBOSE_OUTER_ITERATIONS, true)


--############################################### Initialize and Execute Solver
chiSolverInitialize(phys1)

chiLBTSSetProperty(phys1, "TIMESTEP", 1e-1)
chiLBTSSetProperty(phys1, "VERBOSITY_LEVEL", 0)
chiLBTSSetProperty(phys1, "TIMESTEP_METHOD", "CRANK_NICHOLSON")

phys1name = chiSolverGetName(phys1);

--for k=1,2 do
--    --chiLBTSSetProperty(phys1, "INHIBIT_ADVANCE", true)
--    chiSolverStep(phys1)
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
    chiSolverStep(phys1)
    FRf = chiLBSComputeFissionRate(phys1,"NEW")
    FRi = chiLBSComputeFissionRate(phys1,"OLD")
    dt = chiLBTSGetProperty(phys1, "TIMESTEP")
    time = chiLBTSGetProperty(phys1, "TIME")
    period = dt/math.log(FRf/FRi)
    chiLog(LOG_0, string.format("%s %4d time=%10.3g dt=%10.4g period=%10.3g FR=%10.3e",
            phys1name,k,time,dt,period,FRf))
end
