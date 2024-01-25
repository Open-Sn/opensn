-- 2D Transport test with localized material source Adjoint generation
-- SDM: PWLD
-- Test: None
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
    nodes[i + 1] =  i * ds
end
meshgen = chi_mesh.OrthogonalMeshGenerator.Create({ node_sets = { nodes, nodes } })
chi_mesh.MeshGenerator.Execute(meshgen)

--############################################### Set Material IDs
vol0 = chi_mesh.RPPLogicalVolume.Create({ infx = true,
                                          infy = true,
                                          infz = true })
chiVolumeMesherSetProperty(MATID_FROMLOGICAL, vol0, 0)

vol1 = chi_mesh.RPPLogicalVolume.Create({ ymin = 0.0, ymax = 0.8 * L,
                                          infx = true,
                                          infz = true })
chiVolumeMesherSetProperty(MATID_FROMLOGICAL, vol1, 1)

vol0b = chi_mesh.RPPLogicalVolume.Create({ xmin = -0.166666 + 2.5, xmax = 0.166666 + 2.5,
                                           infy = true,
                                           infz = true })
chiVolumeMesherSetProperty(MATID_FROMLOGICAL, vol0b, 0)

vol2 = chi_mesh.RPPLogicalVolume.Create({ xmin = -0.166666 + 2.5, xmax = 0.166666 + 2.5,
                                          ymin = 0.0, ymax = 2 * 0.166666,
                                          infz = true })
chiVolumeMesherSetProperty(MATID_FROMLOGICAL, vol2, 2)

vol1b = chi_mesh.RPPLogicalVolume.Create({ xmin = -1 + 2.5, xmax = 1 + 2.5,
                                           ymin = 0.9 * L, ymax = L,
                                           infz = true })
chiVolumeMesherSetProperty(MATID_FROMLOGICAL, vol1b, 1)


--############################################### Add materials
num_groups = 1

materials = {}
materials[1] = chiPhysicsAddMaterial("Test Material");
materials[2] = chiPhysicsAddMaterial("Test Material2");
materials[3] = chiPhysicsAddMaterial("Test Material3");

-- Cross sections
chiPhysicsMaterialAddProperty(materials[1], TRANSPORT_XSECTIONS)
chiPhysicsMaterialSetProperty(materials[1], TRANSPORT_XSECTIONS,
                              SIMPLEXS1, 1, 0.01, 0.01)

chiPhysicsMaterialAddProperty(materials[2], TRANSPORT_XSECTIONS)
chiPhysicsMaterialSetProperty(materials[2], TRANSPORT_XSECTIONS,
                              SIMPLEXS1, 1, 0.1 * 20, 0.8)

chiPhysicsMaterialAddProperty(materials[3], TRANSPORT_XSECTIONS)
chiPhysicsMaterialSetProperty(materials[3], TRANSPORT_XSECTIONS,
                              SIMPLEXS1, num_groups, 0.3 * 20, 0.0)

--############################################### Setup Physics
pquad0 = chiCreateProductQuadrature(GAUSS_LEGENDRE_CHEBYSHEV, 48, 6)
chiOptimizeAngularQuadratureForPolarSymmetry(pquad0, 4.0 * math.pi)

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

--############################################### Create QoI region
qoi_vol = chi_mesh.RPPLogicalVolume.Create({ xmin = 0.5, xmax = 0.8333,
                                             ymin = 4.16666, ymax = 4.33333,
                                             infz = true })

chiAdjointSolverAddResponseFunction(phys, "QoI", qoi_vol)
chiSolverSetBasicOption(phys, "REFERENCE_RF", "QoI")

--############################################### Initialize and Execute Solver
ss_solver = lbs.SteadyStateSolver.Create({ lbs_solver_handle = phys })

chiSolverInitialize(ss_solver)
chiSolverExecute(ss_solver)

chiLBSWriteFluxMoments(phys, "Adjoint2D_1b_adjoint")

--############################################### Exports
if master_export == nil then
    ff_m0 = chiGetFieldFunctionHandleByName("phi_g000_m00")
    ff_m1 = chiGetFieldFunctionHandleByName("phi_g000_m01")
    ff_m2 = chiGetFieldFunctionHandleByName("phi_g000_m02")
    chiExportMultiFieldFunctionToVTK({ ff_m0, ff_m1, ff_m2 }, "ZPhi_LBAdjoint")
end
