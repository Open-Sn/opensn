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
N=2000
L=100.0
xmin = -L/2
dx = L/N
for i=1,(N+1) do
    k=i-1
    nodes[i] = xmin + k*dx
end

meshgen1 = mesh.OrthogonalMeshGenerator.Create({ node_sets = {nodes} })
mesh.MeshGenerator.Execute(meshgen1)

--############################################### Set Material IDs
VolumeMesherSetMatIDToAll(0)

vol0 = mesh.RPPLogicalVolume.Create({infx=true, infy=true, zmin=-L/4, zmax=L/4})
VolumeMesherSetProperty(MATID_FROMLOGICAL,vol0,1)

--############################################### Add materials
materials = {}
materials[1] = PhysicsAddMaterial("Strong fuel");
materials[2] = PhysicsAddMaterial("Weak fuel");

PhysicsMaterialAddProperty(materials[1],TRANSPORT_XSECTIONS)
PhysicsMaterialAddProperty(materials[2],TRANSPORT_XSECTIONS)
PhysicsMaterialAddProperty(materials[1],ISOTROPIC_MG_SOURCE)
PhysicsMaterialAddProperty(materials[2],ISOTROPIC_MG_SOURCE)

-- Define microscopic cross sections
xs_strong_fuel_micro = chiPhysicsTransportXSCreate()
chiPhysicsTransportXSSet(xs_strong_fuel_micro, CHI_XSFILE, "tests/transport_transient/xs_inf_k1_6_1g.cxs")
xs_weak_fuelA_micro = chiPhysicsTransportXSCreate()
chiPhysicsTransportXSSet(xs_weak_fuelA_micro, CHI_XSFILE, "tests/transport_transient/xs_inf_critical_1g.cxs")
xs_weak_fuelB_micro = chiPhysicsTransportXSCreate()
chiPhysicsTransportXSSet(xs_weak_fuelB_micro, CHI_XSFILE, "tests/transport_transient/xs_inf_weak_1g.cxs")

atom_density = 0.056559
xs_strong_fuel = chiPhysicsTransportXSMakeCombined({{xs_strong_fuel_micro, atom_density}}) --critical
xs_weak_fuelA  = chiPhysicsTransportXSMakeCombined({{ xs_weak_fuelA_micro, atom_density}}) --critical
xs_weak_fuelB  = chiPhysicsTransportXSMakeCombined({{xs_weak_fuelB_micro, atom_density}})   --critical

num_groups = 1
PhysicsMaterialSetProperty(materials[1],TRANSPORT_XSECTIONS,
        EXISTING,xs_strong_fuel)
PhysicsMaterialSetProperty(materials[2],TRANSPORT_XSECTIONS,
        EXISTING,xs_weak_fuelA)

src={0.0}
PhysicsMaterialSetProperty(materials[1],ISOTROPIC_MG_SOURCE,FROM_ARRAY,src)
PhysicsMaterialSetProperty(materials[2],ISOTROPIC_MG_SOURCE,FROM_ARRAY,src)

function SwapXS(solver_handle, new_xs)
    PhysicsMaterialSetProperty(materials[2],TRANSPORT_XSECTIONS,
            EXISTING,new_xs)
    chiLBSInitializeMaterials(solver_handle)
end

--############################################### Setup Physics
phys1 = chiLBSCreateTransientSolver()

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
--LBSGroupsetSetIterativeMethod(phys1,cur_gs,KRYLOV_RICHARDSON)
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
--LBSSetProperty(phys1,BOUNDARY_CONDITION,ZMIN,LBSBoundaryTypes.REFLECTING);
--LBSSetProperty(phys1,BOUNDARY_CONDITION,ZMAX,LBSBoundaryTypes.REFLECTING);
--
LBSSetProperty(phys1,DISCRETIZATION_METHOD,PWLD)
LBSSetProperty(phys1,SCATTERING_ORDER,1)

chiLBKESSetProperty(phys1, "MAX_ITERATIONS", 1000)
chiLBKESSetProperty(phys1, "TOLERANCE", 1.0e-8)

LBSSetProperty(phys1, USE_PRECURSORS, true)

--LBSSetProperty(phys1, VERBOSE_INNER_ITERATIONS, false)
LBSSetProperty(phys1, VERBOSE_INNER_ITERATIONS, false)
LBSSetProperty(phys1, VERBOSE_OUTER_ITERATIONS, true)


--############################################### Initialize and Execute Solver
SolverInitialize(phys1)

chiLBTSSetProperty(phys1, "TIMESTEP", 1e-3*100)
chiLBTSSetProperty(phys1, "TIMESTOP", 1.0*100)
chiLBTSSetProperty(phys1, "MAX_TIMESTEPS", -1)
chiLBTSSetProperty(phys1, "VERBOSITY_LEVEL", 0)
chiLBTSSetProperty(phys1, "TIMESTEP_METHOD", "CRANK_NICHOLSON")

phys1name = chiSolverGetName(phys1);
initial_FR = chiLBSComputeFissionRate(phys1,"OLD")


--time = 0.0
--psi_t = psi_0
--for k=1,10 do
--    SolverStep(phys1)
--
--    FRf = chiLBSComputeFissionRate(phys1,"NEW") --time+dt
--    FRi = chiLBSComputeFissionRate(phys1,"OLD") --time
--    dt = chiLBTSGetProperty(phys1, "TIMESTEP")
--    time = chiLBTSGetProperty(phys1, "TIME")
--    new_time = time+dt
--
--    period = dt/math.log(FRf/FRi)
--    chiLog(LOG_0, string.format("%s %4d time=%10.3g dt=%10.4g period=%10.3g FR=%10.3e",
--            phys1name,k,time,dt,period,FRf/initial_FR))
--
--    if (not timestep_rejected) then
--        chiLBTSAdvanceTimeValues()
--    end
--    chiLBTSAdvanceTimeData(phys1)
--end

--time = 0.0
--time_stop = 1.0*100
--k=0
--swapped = false
--timestep_rejected = false
--
--tolA = 10.0
--while (time < time_stop) do
--    SolverStep(phys1)
--    FRf = chiLBSComputeFissionRate(phys1,"NEW")
--    FRi = chiLBSComputeFissionRate(phys1,"OLD")
--    dt = chiLBTSGetProperty(phys1, "TIMESTEP")
--    time = chiLBTSGetProperty(phys1, "TIME")
--    period = dt/math.log(FRf/FRi)
--
--    if (time >= 0.2 and not swapped) then
--        SwapXS(phys1, xs_weak_fuelB)
--        swapped = true
--    end
--
--    if (not timestep_rejected) then
--        chiLBTSAdvanceTimeData(phys1)
--        k = k + 1
--        chiLog(LOG_0, string.format("%s %4d time=%10.3g dt=%10.4g period=%10.3g FR=%10.3e",
--                phys1name,k,time,dt,period,FRf/initial_FR))
--    else
--        timestep_rejected = false
--    end
--end
swapped=false
function MyCallBack()
    FRf = chiLBSComputeFissionRate(phys1,"NEW")
    FRi = chiLBSComputeFissionRate(phys1,"OLD")
    dt = chiLBTSGetProperty(phys1, "TIMESTEP")
    time = chiLBTSGetProperty(phys1, "TIME")
    period = dt/math.log(FRf/FRi)

    if (time >= 0.2 and not swapped) then
        SwapXS(phys1, xs_weak_fuelB)
        swapped = true
    end
    chiLog(LOG_0,string.format("%s time=%10.3g dt=%10.4g period=%10.3g FR=%10.3e",
                phys1name,time,dt,period,FRf/initial_FR))
end
chiLBTSSetProperty(phys1, "CALLBACK", "MyCallBack")
chiSolverExecute(phys1)
