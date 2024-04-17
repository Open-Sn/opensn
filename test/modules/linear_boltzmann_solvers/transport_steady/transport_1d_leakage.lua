-- 1D Transport leakage test
-- Unit angular flux left boundary condition in a pure absorber with unit
-- length and a unit absorption cross section. The analytic solution is:
-- j^+ = \int_{0}^{1} \mu e^{-1/\mu} d\mu = 0.10969

-- Check num_procs
num_procs = 3
if (check_num_procs == nil and number_of_processes ~= num_procs) then
    log.Log(LOG_0ERROR, "Incorrect amount of processors. " ..
        "Expected " .. tostring(num_procs) ..
        ". Pass check_num_procs=false to override if possible.")
    os.exit(false)
end

-- Setup mesh
N = 100
L = 1.0
nodes = {}
for i = 1, (N + 1) do
    k = i - 1
    nodes[i] = (i - 1) * L / N
end

meshgen = mesh.OrthogonalMeshGenerator.Create({ node_sets = { nodes } })
mesh.MeshGenerator.Execute(meshgen)
mesh.SetUniformMaterialID(0)

-- Add materials
num_groups = 1
sigma_t = 1.0

materials = {}
materials[1] = mat.AddMaterial("Test Material");
mat.AddProperty(materials[1], TRANSPORT_XSECTIONS)

mat.SetProperty(materials[1], TRANSPORT_XSECTIONS, MODEL1G, sigma_t, 0.0)

-- Setup Physics
pquad = aquad.CreateProductQuadrature(GAUSS_LEGENDRE, 128)
lbs_block = {
    num_groups = num_groups,
    groupsets = {
        {
            groups_from_to = { 0, num_groups - 1 },
            angular_quadrature_handle = pquad,
            angle_aggregation_num_subsets = 1,
            groupset_num_subsets = 1,
            inner_linear_method = "gmres",
            l_abs_tol = 1.0e-6,
            l_max_its = 300,
            gmres_restart_interval = 100,
        }
    }
}

bsrc = {}
for g = 1, num_groups do
    bsrc[g] = 0.0
end
bsrc[1] = 1.0

lbs_options = {
    boundary_conditions = {
        {
            name = "zmin",
            type = "isotropic",
            group_strength = bsrc
        }
    },
    scattering_order = 0,
    save_angular_flux =  true
}

phys = lbs.DiscreteOrdinatesSolver.Create(lbs_block)
lbs.SetOptions(phys, lbs_options)

ss_solver = lbs.SteadyStateSolver.Create({ lbs_solver_handle = phys })

-- Solve the problem
solver.Initialize(ss_solver)
solver.Execute(ss_solver)

-- Compute the leakage
leakage = lbs.ComputeLeakage(phys)
for k, v in pairs(leakage) do
    log.Log(LOG_0, string.format("%s=%.5e", k, v[1]))
end
