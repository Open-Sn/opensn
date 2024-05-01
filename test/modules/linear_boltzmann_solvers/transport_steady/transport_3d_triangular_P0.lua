-- 3D LinearBSolver Vacuum with isotropic incident, triangular set
num_procs = 8

if (check_num_procs==nil and number_of_processes ~= num_procs) then
    log.Log(LOG_0ERROR,"Incorrect amount of processors. " ..
            "Expected "..tostring(num_procs)..
            ". Pass check_num_procs=false to override if possible.")
    os.exit(false)
end

nodes = {}
N = 65
ds = 10.0/N
for i = 0, N do
    nodes[i+1] = i*ds
end
meshgen1 = mesh.OrthogonalMeshGenerator.Create({ node_sets = {nodes,nodes,nodes} })
mesh.MeshGenerator.Execute(meshgen1)
vol0 = logvol.RPPLogicalVolume.Create({infx=true, infy=true, infz=true})
mesh.SetMaterialIDFromLogicalVolume(vol0,0)
materials = {}
materials[1] = mat.AddMaterial("Test Material");
mat.AddProperty(materials[1], TRANSPORT_XSECTIONS)

num_groups = 1
mat.SetProperty(materials[1], TRANSPORT_XSECTIONS, OPENSN_XSFILE, "simple_scatter.xs")
method = 3
sn = 4
moments = 0
tquad = aquad.CreateTriangleQuadrature(method,sn)

lbs_block =
{
    num_groups = num_groups,
    groupsets =
    {
        {
            groups_from_to = {0,0},
            angular_quadrature_handle = tquad,
            angle_aggregation_type = "single",
            angle_aggregation_num_subsets = 1,
            groupset_num_subsets = 1,
            inner_linear_method = "gmres",
            l_abs_tol = 1.0e-8,
            l_max_its = 300,
            gmres_restart_interval = 10,
        },
    }
}

lbs_options =
{
    boundary_conditions =
    {
        {
            name = "xmin",
            type = "isotropic",
            group_strength={1.0},
        }
    },
    scattering_order = 0,
}

phys1 = lbs.DiscreteOrdinatesSolver.Create(lbs_block)
lbs.SetOptions(phys1, lbs_options)
ss_solver = lbs.SteadyStateSolver.Create({lbs_solver_handle = phys1})
solver.Initialize(ss_solver)
solver.Execute(ss_solver)
fflist,count = lbs.GetScalarFieldFunctionList(phys1)
pp1 = post.CellVolumeIntegralPostProcessor.Create
({
  name="max-grp0",
  field_function = fflist[1],
  compute_volume_average = true,
  print_numeric_format = "scientific"
})
post.Execute({ pp1 })
if (master_export == nil) then
    fieldfunc.ExportToVTKMulti(fflist,"ZPhi3D")
end


