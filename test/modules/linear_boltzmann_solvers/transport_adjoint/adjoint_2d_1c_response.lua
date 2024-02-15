-- 2D Transport test with localized material source Adjoint response
-- SDM: PWLD
-- Test: Inner-product=1.38405e-05
num_procs = 4

--############################################### Check num_procs
if (check_num_procs == nil and number_of_processes ~= num_procs) then
    Log(LOG_0ERROR, "Incorrect amount of processors. " ..
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
mesh.SetUniformMaterialID(0)

vol1a = mesh.RPPLogicalVolume.Create(
    {
        infx = true,
        ymin = 0.0, ymax = 0.8 * L,
        infz = true
    }
)

mesh.SetProperty(MATID_FROMLOGICAL, vol1a, 1)

vol0 = mesh.RPPLogicalVolume.Create(
    {
        xmin = 2.5 - 0.166666, xmax = 2.5 + 0.166666,
        infy = true,
        infz = true
    }
)
mesh.SetProperty(MATID_FROMLOGICAL, vol0, 0)

vol2 = mesh.RPPLogicalVolume.Create(
    {
        xmin = 2.5 - 0.166666, xmax = 2.5 + 0.166666,
        ymin = 0.0, ymax = 2 * 0.166666,
        infz = true
    }
)
mesh.SetProperty(MATID_FROMLOGICAL, vol2, 2)

vol1b = mesh.RPPLogicalVolume.Create(
    {
        xmin = -1 + 2.5, xmax = 1 + 2.5,
        ymin = 0.9 * L, ymax = L,
        infz = true
    }
)
mesh.SetProperty(MATID_FROMLOGICAL, vol1b, 1)

--############################################### Add materials
num_groups = 1

materials = {}
materials[1] = PhysicsAddMaterial("Test Material");
materials[2] = PhysicsAddMaterial("Test Material2");
materials[3] = PhysicsAddMaterial("Test Material3");

-- Cross sections
PhysicsMaterialAddProperty(materials[1], TRANSPORT_XSECTIONS)
PhysicsMaterialSetProperty(materials[1], TRANSPORT_XSECTIONS,
                           SIMPLEXS1, 1, 0.01, 0.01)

PhysicsMaterialAddProperty(materials[2], TRANSPORT_XSECTIONS)
PhysicsMaterialSetProperty(materials[2], TRANSPORT_XSECTIONS,
                           SIMPLEXS1, 1, 0.1 * 20, 0.8)

PhysicsMaterialAddProperty(materials[3], TRANSPORT_XSECTIONS)
PhysicsMaterialSetProperty(materials[3], TRANSPORT_XSECTIONS,
                           SIMPLEXS1, num_groups, 0.3 * 20, 0.0)

-- Sources
src = {}
for g = 1, num_groups do
    if g == 1 then src[g] = 3.0 else src[g] = 0.0 end
end
PhysicsMaterialAddProperty(materials[3], ISOTROPIC_MG_SOURCE)
PhysicsMaterialSetProperty(materials[3], ISOTROPIC_MG_SOURCE, FROM_ARRAY, src)

--############################################### Setup Physics
pquad0 = CreateProductQuadrature(GAUSS_LEGENDRE_CHEBYSHEV, 48, 6)
OptimizeAngularQuadratureForPolarSymmetry(pquad0, 4.0 * math.pi)

lbs_block = {
    num_groups = num_groups,
    groupsets = {
        {
            groups_from_to = { 0, num_groups - 1 },
            angular_quadrature_handle = pquad0,
            inner_linear_method = "gmres",
            l_abs_tol = 1.0e-6,
            l_max_its = 500,
            gmres_restart_interval = 100,
        },
    }
}

lbs_options = { scattering_order = 1 }

phys = lbs.DiscreteOrdinatesAdjointSolver.Create(lbs_block)
lbs.SetOptions(phys, lbs_options)

--############################################### Initialize and Compute Response
ss_solver = lbs.SteadyStateSolver.Create({ lbs_solver_handle = phys })

SolverInitialize(ss_solver)
LBSReadFluxMoments(phys, "Adjoint2D_1b_adjoint")

value = AdjointSolverComputeInnerProduct(phys)
Log(LOG_0, string.format("Inner-product=%.5e", value))

--############################################### Exports
if master_export == nil then
    ff_m0 = GetFieldFunctionHandleByName("phi_g000_m00")
    ff_m1 = GetFieldFunctionHandleByName("phi_g000_m01")
    ff_m2 = GetFieldFunctionHandleByName("phi_g000_m02")

    ExportMultiFieldFunctionToVTK({ ff_m0, ff_m1, ff_m2 }, "ZPhi_LBAdjointResponse")
end

--############################################### Cleanup
MPIBarrier()
if (location_id == 0) then
    os.execute("rm Adjoint2D_1b_adjoint*.data")
end
