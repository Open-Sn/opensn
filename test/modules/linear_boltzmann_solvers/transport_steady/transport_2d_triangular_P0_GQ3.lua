-- 2D LinearBSolver Vacuum with isotropic incident, triangular set, Henyey-GS Forward peaked cross-sections
-- SDM: PWLD
-- Test: 
num_procs = 1

--############################################### Check num_procs
if (check_num_procs==nil and number_of_processes ~= num_procs) then
    log.Log(LOG_0ERROR,"Incorrect amount of processors. " ..
            "Expected "..tostring(num_procs)..
            ". Pass check_num_procs=false to override if possible.")
    os.exit(false)
end

--############################################### Setup mesh
nodes = {}
N = 65
ds = 1.0 / N
for i = 0, N do
    nodes[i + 1] = i * ds
end
meshgen1 = mesh.OrthogonalMeshGenerator.Create({ node_sets = { nodes, nodes } })
mesh.MeshGenerator.Execute(meshgen1)

--############################################### Set Material IDs
vol0 = logvol.RPPLogicalVolume.Create({ infx = true, infy = true, infz = true })
mesh.SetMaterialIDFromLogicalVolume(vol0, 0)
--############################################### Add materials
materials = {}
materials[1] = mat.AddMaterial("Test Material");

mat.AddProperty(materials[1], TRANSPORT_XSECTIONS)

mat.AddProperty(materials[1], ISOTROPIC_MG_SOURCE)

num_groups = 1
mat.SetProperty(materials[1], TRANSPORT_XSECTIONS, OPENSN_XSFILE, "xs_Henyey_GS_p4.xs")

src = {}
for g = 1, num_groups do
    src[g] = 0.0
end
--src[1] = 1.0
mat.SetProperty(materials[1], ISOTROPIC_MG_SOURCE, FROM_ARRAY, src)

method = 3
scattering_order = 0
sn = 8
moments = 0

pquad = aquad.CreateTriangleQuadrature(method, sn,scattering_order)
aquad.OptimizeForPolarSymmetry(pquad,4.0*math.pi)

--========== Groupset def

lbs_block = {
    num_groups = num_groups,
    groupsets = {
        {
            groups_from_to = { 0, 0 },
            angular_quadrature_handle = pquad,
            angle_aggregation_type = "single",
            angle_aggregation_num_subsets = 1,
            groupset_num_subsets = 1,
            inner_linear_method = "gmres",
            l_abs_tol = 1.0e-8,
            l_max_its = 300,
            gmres_restart_interval = 50,
        },
    }
}

lbs_options = {
    boundary_conditions = {
        {
            name = "xmin",
            type = "isotropic",
            group_strength = { 1.0 },
        }
    },
    scattering_order = 0,
}

phys1 = lbs.DiscreteOrdinatesSolver.Create(lbs_block)
lbs.SetOptions(phys1, lbs_options)

--========== Solvers
ss_solver = lbs.SteadyStateSolver.Create({ lbs_solver_handle = phys1 })

solver.Initialize(ss_solver)
solver.Execute(ss_solver)

--############################################### Setup Output
ss_solver = lbs.SteadyStateSolver.Create({ lbs_solver_handle = phys1 })

solver.Initialize(ss_solver)
solver.Execute(ss_solver)

--############################################### Get field functions
fflist, count = lbs.GetScalarFieldFunctionList(phys1)

--############################################### Slice plot
--slice2 = fieldfunc.FFInterpolationCreate(SLICE)
--fieldfunc.SetProperty(slice2, SLICE_POINT, 0.0, 0.0, 0.025)
--fieldfunc.SetProperty(slice2, ADD_FIELDFUNCTION, fflist[1])

--fieldfunc.Initialize(slice2)
--fieldfunc.Execute(slice2)

--############################################### Volume integrations
--ffi1 = fieldfunc.FFInterpolationCreate(VOLUME)
--curffi = ffi1
--fieldfunc.SetProperty(curffi, OPERATION, OP_MAX)
--fieldfunc.SetProperty(curffi, LOGICAL_VOLUME, vol0)
--fieldfunc.SetProperty(curffi, ADD_FIELDFUNCTION, fflist[1])

--fieldfunc.Initialize(curffi)
--fieldfunc.Execute(curffi)
--maxval = fieldfunc.GetValue(curffi)

--log.Log(LOG_0, string.format("Max-value1=%.5f", maxval))

--############################################### Exports
if master_export == nil then
    fieldfunc.ExportToVTKMulti(fflist,"ZPhi")
end

