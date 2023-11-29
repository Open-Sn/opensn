-- 2D Transport test with point source Multigroup Adjoint generation
-- SDM: PWLD
-- Test: None
num_procs = 4

--############################################### Check num_procs
if (check_num_procs == nil and number_of_processes ~= num_procs) then
    chiLog(LOG_0ERROR, "Incorrect amount of processors. " ..
        "Expected " .. tostring(num_procs) ..
        ". Pass check_num_procs=false to override if possible.")
    os.exit(false)
end

--############################################### Setup mesh
N = 60
L = 5.0
ds = L / N

nodes = {}
for i = 0, N do
    nodes[i + 1] = i * ds
end
meshgen = mesh.OrthogonalMeshGenerator.Create({ node_sets = { nodes, nodes } })
mesh.MeshGenerator.Execute(meshgen)

--############################################### Set Material IDs
VolumeMesherSetMatIDToAll(0)

vol1a = mesh.RPPLogicalVolume.Create(
    {
        infx = true,
        ymin = 0.0, ymax = 0.8 * L,
        infz = true
    }
)

VolumeMesherSetProperty(MATID_FROMLOGICAL, vol1a, 1)

vol0 = mesh.RPPLogicalVolume.Create(
    {
        xmin = 2.5 - 0.166666, xmax = 2.5 + 0.166666,
        infy = true,
        infz = true
    }
)
VolumeMesherSetProperty(MATID_FROMLOGICAL, vol0, 0)

vol1b = mesh.RPPLogicalVolume.Create(
    {
        xmin = -1 + 2.5, xmax = 1 + 2.5,
        ymin = 0.9 * L, ymax = L,
        infz = true
    }
)
VolumeMesherSetProperty(MATID_FROMLOGICAL, vol1b, 1)

--############################################### Add materials
num_groups = 10

materials = {}
materials[1] = PhysicsAddMaterial("Test Material1");
materials[2] = PhysicsAddMaterial("Test Material2");

-- Cross sections
PhysicsMaterialAddProperty(materials[1], TRANSPORT_XSECTIONS)
PhysicsMaterialSetProperty(materials[1], TRANSPORT_XSECTIONS,
                           SIMPLEXS1, num_groups, 0.01, 0.01)

PhysicsMaterialAddProperty(materials[2], TRANSPORT_XSECTIONS)
PhysicsMaterialSetProperty(materials[2], TRANSPORT_XSECTIONS,
                           SIMPLEXS1, num_groups, 0.1 * 20, 0.8)

-- Distributed sources
qoi_vol = mesh.RPPLogicalVolume.Create(
    {
        xmin = 0.5, xmax = 0.8333,
        ymin = 4.16666, ymax = 4.33333,
        infz = true
    }
)

function ResponseFunction(xyz, mat_id)
    response = {}
    for g = 1, num_groups do
        if g == 6 then response[g] = 1.0 else response[g] = 0.0 end
    end
    return response
end
response_func = opensn.LuaSpatialMaterialFunction.Create({ lua_function_name = "ResponseFunction" })

dsrc = lbs.DistributedSource.Create(
    {
        logical_volume_handle = qoi_vol,
        function_handle = response_func
    }
)

--############################################### Setup Physics
pquad = CreateProductQuadrature(GAUSS_LEGENDRE_CHEBYSHEV, 12, 2)
chiOptimizeAngularQuadratureForPolarSymmetry(pquad, 4.0 * math.pi)

lbs_block = {
    num_groups = num_groups,
    groupsets = {
        {
            groups_from_to = { 0, num_groups - 1 },
            angular_quadrature_handle = pquad,
            inner_linear_method = "gmres",
            l_abs_tol = 1.0e-6,
            l_max_its = 500,
            gmres_restart_interval = 100,
        },
    }
}

lbs_options = {
    scattering_order = 1,
    adjoint = true,
    distributed_sources = { dsrc }
}

phys = lbs.DiscreteOrdinatesSolver.Create(lbs_block)
lbs.SetOptions(phys, lbs_options)

--############################################### Initialize and Execute Solver
ss_solver = lbs.SteadyStateSolver.Create({ lbs_solver_handle = phys })

SolverInitialize(ss_solver)
chiSolverExecute(ss_solver)

chiLBSWriteFluxMoments(phys, "Adjoint2D_3b_adjoint")

--############################################### Exports
if master_export == nil then
    ff_m0 = chiGetFieldFunctionHandleByName("phi_g000_m00")
    ff_m1 = chiGetFieldFunctionHandleByName("phi_g000_m01")
    ff_m2 = chiGetFieldFunctionHandleByName("phi_g000_m02")
    chiExportMultiFieldFunctionToVTK({ ff_m0, ff_m1, ff_m2 }, "ZPhi_LBAdjoint")
end
