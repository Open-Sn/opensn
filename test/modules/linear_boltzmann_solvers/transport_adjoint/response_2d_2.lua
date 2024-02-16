-- 2D Transport test with point source FWD
-- SDM: PWLD
-- Test: QoI Value=2.90386e-05
--       Inner Product=2.90543e-05
num_procs = 4

-- Check num_procs
if (check_num_procs == nil and number_of_processes ~= num_procs) then
    Log(LOG_0ERROR, "Incorrect amount of processors. " ..
            "Expected " .. tostring(num_procs) ..
            ". Pass check_num_procs=false to override if possible.")
    os.exit(false)
end

-- Setup mesh
N = 60
L = 5.0
ds = L / N

nodes = {}
for i = 0, N do
    nodes[i + 1] = i * ds
end
meshgen = mesh.OrthogonalMeshGenerator.Create({ node_sets = { nodes, nodes } })
mesh.MeshGenerator.Execute(meshgen)

-- Set Material IDs
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

vol1b = mesh.RPPLogicalVolume.Create(
        {
            xmin = -1 + 2.5, xmax = 1 + 2.5,
            ymin = 0.9 * L, ymax = L,
            infz = true
        }
)
mesh.SetProperty(MATID_FROMLOGICAL, vol1b, 1)

-- Add materials
num_groups = 1

materials = {}
materials[1] = PhysicsAddMaterial("Test Material1");
materials[2] = PhysicsAddMaterial("Test Material2");

-- Add cross sections
PhysicsMaterialAddProperty(materials[1], TRANSPORT_XSECTIONS)
PhysicsMaterialSetProperty(materials[1], TRANSPORT_XSECTIONS,
        SIMPLEXS1, num_groups, 0.01, 0.01)

PhysicsMaterialAddProperty(materials[2], TRANSPORT_XSECTIONS)
PhysicsMaterialSetProperty(materials[2], TRANSPORT_XSECTIONS,
        SIMPLEXS1, num_groups, 0.1 * 20, 0.8)

-- Add sources
src = {}
for g = 1, num_groups do
    if g == 1 then
        src[g] = 1.0
    else
        src[g] = 0.0
    end
end

loc = { 1.25 - 0.5 * ds, 1.5 * ds, 0.0 }
pt_src = lbs.PointSource.Create({ location = loc, strength = src })

-- Setup physics
pquad = CreateProductQuadrature(GAUSS_LEGENDRE_CHEBYSHEV, 48, 6)
OptimizeAngularQuadratureForPolarSymmetry(pquad, 4.0 * math.pi)

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
        }
    },
    options = {
        scattering_order = 0,
        point_sources = { pt_src }
    }
}
phys = lbs.DiscreteOrdinatesSolver.Create(lbs_block)

-- Forward solve
ss_solver = lbs.SteadyStateSolver.Create({ lbs_solver_handle = phys })

SolverInitialize(ss_solver)
SolverExecute(ss_solver)

-- Get field functions
ff_m0 = GetFieldFunctionHandleByName("phi_g000_m00")
ff_m1 = GetFieldFunctionHandleByName("phi_g000_m01")
ff_m2 = GetFieldFunctionHandleByName("phi_g000_m02")

-- Define QoI region
qoi_vol = mesh.RPPLogicalVolume.Create(
        {
            xmin = 0.5, xmax = 0.8333,
            ymin = 4.16666, ymax = 4.33333,
            infz = true
        }
)

-- Compute QoI
ffi = FFInterpolationCreate(VOLUME)
FFInterpolationSetProperty(ffi, OPERATION, OP_SUM)
FFInterpolationSetProperty(ffi, LOGICAL_VOLUME, qoi_vol)
FFInterpolationSetProperty(ffi, ADD_FIELDFUNCTION, ff_m0)

FFInterpolationInitialize(ffi)
FFInterpolationExecute(ffi)
fwd_qoi = FFInterpolationGetValue(ffi)

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
SolverExecute(ss_solver)
LBSWriteFluxMoments(phys, "adjoint_2d_2")

-- Create response evaluator
buffers = { { name = "buff", file_prefixes = { flux_moments = "adjoint_2d_2" } } }
pt_sources = { pt_src }
response_options = {
    lbs_solver_handle = phys,
    options = {
        buffers = buffers,
        sources = { point = pt_sources }
    }
}
evaluator = lbs.ResponseEvaluator.Create(response_options)

-- Evaluate response
adj_qoi = lbs.EvaluateResponse(evaluator, "buff")

-- Print results
Log(LOG_0, string.format("QoI Value=%.5e", fwd_qoi))
Log(LOG_0, string.format("Inner Product=%.5e", adj_qoi))

-- Cleanup
MPIBarrier()
if (location_id == 0) then
    os.execute("rm adjoint_2d_2*")
end
