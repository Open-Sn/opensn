-- 2D Transport test with point source FWD
-- SDM: PWLD
-- Test: QOI-value=2.90386e-05
num_procs = 4





--############################################### Check num_procs
if (check_num_procs==nil and chi_number_of_processes ~= num_procs) then
    chiLog(LOG_0ERROR,"Incorrect amount of processors. " ..
                      "Expected "..tostring(num_procs)..
                      ". Pass check_num_procs=false to override if possible.")
    os.exit(false)
end

--############################################### Setup mesh
nodes={}
N=60
L=5.0
ds=L/N
xmin=0.0
for i=0,N do
    nodes[i+1] = xmin + i*ds
end
meshgen1 = chi_mesh.OrthogonalMeshGenerator.Create({ node_sets = {nodes,nodes} })
chi_mesh.MeshGenerator.Execute(meshgen1)

----############################################### Set Material IDs
NewRPP = chi_mesh.RPPLogicalVolume.Create
vol0 = NewRPP({infx=true, infy=true, infz=true})
vol1 = NewRPP({ymin=0.0,ymax=0.8*L,infx=true,infz=true})
chiVolumeMesherSetProperty(MATID_FROMLOGICAL,vol0,0)
chiVolumeMesherSetProperty(MATID_FROMLOGICAL,vol1,1)



----############################################### Set Material IDs
vol0b = NewRPP({xmin=-0.166666+2.5,xmax=0.166666+2.5,infy=true,infz=true})
vol1b = NewRPP({xmin=-1+2.5,xmax=1+2.5,ymin=0.9*L,ymax=L,infz=true})
chiVolumeMesherSetProperty(MATID_FROMLOGICAL,vol0b,0)
chiVolumeMesherSetProperty(MATID_FROMLOGICAL,vol1b,1)


--############################################### Add materials
materials = {}
materials[1] = chiPhysicsAddMaterial("Test Material");
materials[2] = chiPhysicsAddMaterial("Test Material2");

chiPhysicsMaterialAddProperty(materials[1],TRANSPORT_XSECTIONS)
chiPhysicsMaterialAddProperty(materials[2],TRANSPORT_XSECTIONS)

chiPhysicsMaterialAddProperty(materials[1],ISOTROPIC_MG_SOURCE)
chiPhysicsMaterialAddProperty(materials[2],ISOTROPIC_MG_SOURCE)


num_groups = 1
chiPhysicsMaterialSetProperty(materials[1],
                              TRANSPORT_XSECTIONS,
                              SIMPLEXS1,1,0.01,0.01)
chiPhysicsMaterialSetProperty(materials[2],
                              TRANSPORT_XSECTIONS,
                              SIMPLEXS1,1,0.1*20,0.8)

src={}
for g=1,num_groups do
    src[g] = 0.0
end
src[1] = 0.0
chiPhysicsMaterialSetProperty(materials[1],ISOTROPIC_MG_SOURCE,FROM_ARRAY,src)
src[1] = 0.0
chiPhysicsMaterialSetProperty(materials[2],ISOTROPIC_MG_SOURCE,FROM_ARRAY,src)
src[1] = 1.0

--############################################### Setup Physics
pquad0 = chiCreateProductQuadrature(GAUSS_LEGENDRE_CHEBYSHEV,48, 6)
chiOptimizeAngularQuadratureForPolarSymmetry(pquad0, 4.0*math.pi)

lbs_block =
{
    num_groups = num_groups,
    groupsets =
    {
        {
            groups_from_to = {0, num_groups-1},
            angular_quadrature_handle = pquad0,
            inner_linear_method = "gmres",
            l_abs_tol = 1.0e-6,
            l_max_its = 500,
            gmres_restart_interval = 100,
        },
    }
}

lbs_options =
{
    scattering_order = 1,
}

phys1 = lbs.DiscreteOrdinatesSolver.Create(lbs_block)
lbs.SetOptions(phys1, lbs_options)

--############################################### Add point source
chiLBSAddPointSource(phys1, 1.25 - 0.5*ds, 1.5*ds, 0.0, src)

--############################################### Initialize and Execute Solver
ss_solver = lbs.SteadyStateSolver.Create({lbs_solver_handle = phys1})

chiSolverInitialize(ss_solver)
chiSolverExecute(ss_solver)

--############################################### Create QOIs
tvol0 = NewRPP({xmin=2.3333,xmax=2.6666,ymin=4.16666,ymax=4.33333,infz=true})
tvol1 = NewRPP({xmin=0.5   ,xmax=0.8333,ymin=4.16666,ymax=4.33333,infz=true})

--############################################### Get field functions
ff_m0 = chiGetFieldFunctionHandleByName("phi_g000_m00")
ff_m1 = chiGetFieldFunctionHandleByName("phi_g000_m01")
ff_m2 = chiGetFieldFunctionHandleByName("phi_g000_m02")


--############################################### Slice plot

--############################################### Volume integrations
ffi1 = chiFFInterpolationCreate(VOLUME)
curffi = ffi1
chiFFInterpolationSetProperty(curffi,OPERATION,OP_SUM)
chiFFInterpolationSetProperty(curffi,LOGICAL_VOLUME,tvol1)
chiFFInterpolationSetProperty(curffi,ADD_FIELDFUNCTION,ff_m0)

chiFFInterpolationInitialize(curffi)
chiFFInterpolationExecute(curffi)
QOI_value = chiFFInterpolationGetValue(curffi)

chiLog(LOG_0,string.format("QOI-value=%.5e", QOI_value))

--############################################### Exports
if master_export == nil then
    chiExportMultiFieldFunctionToVTK({ff_m0, ff_m1, ff_m2},"ZPhi_LBS")
end
