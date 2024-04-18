-- 2D Transport test with localized material source
-- SDM: PWLD
-- Test: QoI Value=1.38399e-05
--       Inner Product=1.38405e-05
num_procs = 4

-- Check num_procs
if (check_num_procs == nil and number_of_processes ~= num_procs) then
    log.Log(LOG_0ERROR, "Incorrect amount of processors. " ..
            "Expected " .. tostring(num_procs) ..
            ". Pass check_num_procs=false to override if possible.")
    os.exit(false)
end

-- Create mesh
N = 60
L = 5.0
ds = L / N

nodes = {}
for i = 0, N do
    nodes[i + 1] = i * ds
end

meshgen = mesh.OrthogonalMeshGenerator.Create({ node_sets = { nodes, nodes } })
mesh.MeshGenerator.Execute(meshgen)

-- Set material IDs
mesh.SetUniformMaterialID(0)

vol1a = logvol.RPPLogicalVolume.Create(
        {
            infx = true,
            ymin = 0.0, ymax = 0.8 * L,
            infz = true
        }
)

mesh.SetMaterialIDFromLogicalVolume(vol1a, 1)

vol0 = logvol.RPPLogicalVolume.Create(
        {
            xmin = 2.5 - 0.166666, xmax = 2.5 + 0.166666,
            infy = true,
            infz = true
        }
)
mesh.SetMaterialIDFromLogicalVolume(vol0, 0)

vol2 = logvol.RPPLogicalVolume.Create(
        {
            xmin = 2.5 - 0.166666, xmax = 2.5 + 0.166666,
            ymin = 0.0, ymax = 2 * 0.166666,
            infz = true
        }
)
mesh.SetMaterialIDFromLogicalVolume(vol2, 2)

vol1b = logvol.RPPLogicalVolume.Create(
        {
            xmin = -1 + 2.5, xmax = 1 + 2.5,
            ymin = 0.9 * L, ymax = L,
            infz = true
        }
)
mesh.SetMaterialIDFromLogicalVolume(vol1b, 1)

-- Create materials
materials = {}
materials[1] = mat.AddMaterial("Test Material1");
materials[2] = mat.AddMaterial("Test Material2");
materials[3] = mat.AddMaterial("Test Material3");

-- Add cross sections to materials
num_groups = 1
mat.AddProperty(materials[1], TRANSPORT_XSECTIONS)
mat.SetProperty(materials[1], TRANSPORT_XSECTIONS, SIMPLE_ONE_GROUP, 0.01, 0.01)

mat.AddProperty(materials[2], TRANSPORT_XSECTIONS)
mat.SetProperty(materials[2], TRANSPORT_XSECTIONS, SIMPLE_ONE_GROUP, 0.1 * 20, 0.8)

mat.AddProperty(materials[3], TRANSPORT_XSECTIONS)
mat.SetProperty(materials[3], TRANSPORT_XSECTIONS, SIMPLE_ONE_GROUP, 0.3 * 20, 0.0)

-- Create sources
src = {}
for g = 1, num_groups do
    if g == 1 then
        src[g] = 3.0
    else
        src[g] = 0.0
    end
end
mat.AddProperty(materials[3], ISOTROPIC_MG_SOURCE)
mat.SetProperty(materials[3], ISOTROPIC_MG_SOURCE, FROM_ARRAY, src)

-- Setup physics
pquad = aquad.CreateProductQuadrature(GAUSS_LEGENDRE_CHEBYSHEV, 48, 6)
aquad.OptimizeForPolarSymmetry(pquad, 4.0 * math.pi)

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
    },
    options = { scattering_order = 0 }
}
phys = lbs.DiscreteOrdinatesSolver.Create(lbs_block)

-- Forward solve
ss_solver = lbs.SteadyStateSolver.Create({ lbs_solver_handle = phys })

solver.Initialize(ss_solver)
solver.Execute(ss_solver)

-- Get field functions
ff_m0 = fieldfunc.GetHandleByName("phi_g000_m00")
ff_m1 = fieldfunc.GetHandleByName("phi_g000_m01")
ff_m2 = fieldfunc.GetHandleByName("phi_g000_m02")

-- Define QoI region
qoi_vol = logvol.RPPLogicalVolume.Create(
        {
            xmin = 0.5, xmax = 0.8333,
            ymin = 4.16666, ymax = 4.33333,
            infz = true
        }
)

-- Compute QoI
ffi = fieldfunc.FFInterpolationCreate(VOLUME)
fieldfunc.SetProperty(ffi, OPERATION, OP_SUM)
fieldfunc.SetProperty(ffi, LOGICAL_VOLUME, qoi_vol)
fieldfunc.SetProperty(ffi, ADD_FIELDFUNCTION, ff_m0)

fieldfunc.Initialize(ffi)
fieldfunc.Execute(ffi)
fwd_qoi = fieldfunc.GetValue(ffi)

-- Create adjoint source
adjoint_source = lbs.DistributedSource.Create(
        { logical_volume_handle = qoi_vol }
)

-- Switch to adjoint mode
adjoint_options = {
    adjoint = true,
    distributed_sources = { adjoint_source }
}
lbs.SetOptions(phys, adjoint_options)

-- Adjoint solve, write results
solver.Execute(ss_solver)
lbs.WriteFluxMoments(phys, "adjoint_2d_1")

-- Create response evaluator
buffers = { { name = "buff", file_prefixes = { flux_moments = "adjoint_2d_1" } } }
mat_sources = { { material_id = 2, strength = src } }
response_options = {
    lbs_solver_handle = phys,
    options = {
        buffers = buffers,
        sources = { material = mat_sources }
    }
}
evaluator = lbs.ResponseEvaluator.Create(response_options)

-- Evaluate response
adj_qoi = lbs.EvaluateResponse(evaluator, "buff")

-- Print results
log.Log(LOG_0, string.format("QoI Value=%.5e", fwd_qoi))
log.Log(LOG_0, string.format("Inner Product=%.5e", adj_qoi))

-- Cleanup
MPIBarrier()
if (location_id == 0) then
    os.execute("rm adjoint_2d_1*")
end
