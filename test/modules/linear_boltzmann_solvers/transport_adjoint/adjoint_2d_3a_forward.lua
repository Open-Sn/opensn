-- 2D Transport test with point source Multigroup FWD
-- SDM: PWLD
-- Test:
--  QOI-value[0]= 1.12687e-06
--  QOI-value[1]= 2.95934e-06
--  QOI-value[2]= 3.92975e-06
--  QOI-value[3]= 4.18474e-06
--  QOI-value[4]= 3.89649e-06
--  QOI-value[5]= 3.30482e-06
--  QOI-value[6]= 1.54506e-06
--  QOI-value[7]= 6.74868e-07
--  QOI-value[8]= 3.06178e-07
--  QOI-value[9]= 2.07284e-07
--  QOI-value[sum]= 2.21354e-05
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

phys = lbs.DiscreteOrdinatesSolver.Create(lbs_block)
lbs.SetOptions(phys, lbs_options)

--############################################### Initialize and Execute Solver
ss_solver = lbs.SteadyStateSolver.Create({ lbs_solver_handle = phys })

chiSolverInitialize(ss_solver)
chiSolverExecute(ss_solver)

--############################################### Create QoI Region
qoi_vol = chi_mesh.RPPLogicalVolume.Create({ xmin = 0.5, xmax = 0.8333,
                                             ymin = 4.16666, ymax = 4.33333,
                                             infz = true })

--############################################### Get field functions
ff_m0 = chiGetFieldFunctionHandleByName("phi_g000_m00")
ff_m1 = chiGetFieldFunctionHandleByName("phi_g000_m01")
ff_m2 = chiGetFieldFunctionHandleByName("phi_g000_m02")


--############################################### Slice plot

--############################################### Volume integrations
sum = 0.0
for g = 0, num_groups - 1 do
    ff = chiGetFieldFunctionHandleByName(
        "phi_g" .. string.format("%03d", g) .. "_m" .. string.format("%02d", 0))
    ffi = chiFFInterpolationCreate(VOLUME)
    chiFFInterpolationSetProperty(ffi, OPERATION, OP_SUM)
    chiFFInterpolationSetProperty(ffi, LOGICAL_VOLUME, qoi_vol)
    chiFFInterpolationSetProperty(ffi, ADD_FIELDFUNCTION, ff)

    chiFFInterpolationInitialize(ffi)
    chiFFInterpolationExecute(ffi)
    value = chiFFInterpolationGetValue(ffi)

    chiLog(LOG_0, string.format("QoI-value[" .. tostring(g) .. "]= %.5e", value))

    sum = sum + value
end
chiLog(LOG_0, string.format("QoI-value[sum]= %.5e", sum))

--############################################### Exports
if master_export == nil then
    chiExportMultiFieldFunctionToVTK({ ff_m0, ff_m1, ff_m2 }, "ZPhi_LBS")
end

--[0]  QOI-vallue[0]= 1.95637e-09
--[0]  QOI-vallue[1]= 5.13773e-09
--[0]  QOI-vallue[2]= 6.82248e-09
--[0]  QOI-vallue[3]= 7.26518e-09
--[0]  QOI-vallue[4]= 6.76479e-09
--[0]  QOI-vallue[5]= 5.73749e-09
--[0]  QOI-vallue[6]= 2.68214e-09
--[0]  QOI-vallue[7]= 1.17145e-09
--[0]  QOI-vallue[8]= 5.31336e-10
--[0]  QOI-vallue[9]= 3.59573e-10
