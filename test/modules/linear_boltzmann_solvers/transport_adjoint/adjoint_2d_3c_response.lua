-- 2D Transport test with point source Multigroup Adjoint Response
-- SDM: PWLD
-- Test: Inner-product=3.30607e-06
num_procs = 4

--############################################### Check num_procs
if (check_num_procs == nil and chi_number_of_processes ~= num_procs) then
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
meshgen = chi_mesh.OrthogonalMeshGenerator.Create({ node_sets = { nodes, nodes } })
chi_mesh.MeshGenerator.Execute(meshgen)

--############################################### Set Material IDs
vol0 = chi_mesh.RPPLogicalVolume.Create({ infx = true,
                                          infy = true,
                                          infz = true })
chiVolumeMesherSetProperty(MATID_FROMLOGICAL, vol0, 0)

vol1 = chi_mesh.RPPLogicalVolume.Create({ infx = true,
                                          ymin = 0.0, ymax = 0.8 * L,
                                          infz = true })
chiVolumeMesherSetProperty(MATID_FROMLOGICAL, vol1, 1)

vol0b = chi_mesh.RPPLogicalVolume.Create({ xmin = -0.166666 + 2.5, xmax = 0.166666 + 2.5,
                                           infy = true,
                                           infz = true })
chiVolumeMesherSetProperty(MATID_FROMLOGICAL, vol0b, 0)

vol1b = chi_mesh.RPPLogicalVolume.Create({ xmin = -1 + 2.5, xmax = 1 + 2.5,
                                           ymin = 0.9 * L, ymax = L,
                                           infz = true })
chiVolumeMesherSetProperty(MATID_FROMLOGICAL, vol1b, 1)


--############################################### Add materials
num_groups = 10

materials = {}
materials[1] = chiPhysicsAddMaterial("Test Material1");
materials[2] = chiPhysicsAddMaterial("Test Material2");

chiPhysicsMaterialAddProperty(materials[1], TRANSPORT_XSECTIONS)
chiPhysicsMaterialSetProperty(materials[1], TRANSPORT_XSECTIONS,
                              SIMPLEXS1, num_groups, 0.01, 0.01)

chiPhysicsMaterialAddProperty(materials[2], TRANSPORT_XSECTIONS)
chiPhysicsMaterialSetProperty(materials[2], TRANSPORT_XSECTIONS,
                              SIMPLEXS1, num_groups, 0.1 * 20, 0.8)

-- Sources
src = {}
for g = 1, num_groups do
    if g == 1 then src[g] = 1.0 else src[g] = 0.0 end
end

loc = { 1.25 - 0.5 * ds, 1.5 * ds, 0.0 }
pt_src = lbs.PointSource.Create({ location = loc, strength = src })

--############################################### Setup Physics
pquad = chiCreateProductQuadrature(GAUSS_LEGENDRE_CHEBYSHEV, 12, 2)
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
    point_sources = {pt_src}
}

phys = lbs.DiscreteOrdinatesAdjointSolver.Create(lbs_block)
lbs.SetOptions(phys, lbs_options)

--############################################### Initialize and Compute Response
ss_solver = lbs.SteadyStateSolver.Create({lbs_solver_handle = phys})

chiSolverInitialize(ss_solver)

chiLBSReadFluxMoments(phys, "Adjoint2D_3b_adjoint")
value = chiAdjointSolverComputeInnerProduct(phys)
chiLog(LOG_0,string.format("Inner-product=%.5e", value))

--############################################### Exports
if master_export == nil then
    ff_m0 = chiGetFieldFunctionHandleByName("phi_g000_m00")
    ff_m1 = chiGetFieldFunctionHandleByName("phi_g000_m01")
    ff_m2 = chiGetFieldFunctionHandleByName("phi_g000_m02")
    chiExportMultiFieldFunctionToVTK({ff_m0, ff_m1, ff_m2},"ZPhi_LBAdjointResponse")
end

--############################################### Cleanup
chiMPIBarrier()
if (chi_location_id == 0) then
    os.execute("rm Adjoint2D_3b_adjoint*.data")
end
