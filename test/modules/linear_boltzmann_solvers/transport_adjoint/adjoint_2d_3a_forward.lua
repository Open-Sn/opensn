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


num_groups = 10
chiPhysicsMaterialSetProperty(materials[1],
                              TRANSPORT_XSECTIONS,
                              SIMPLEXS1,num_groups,0.01,0.01)
chiPhysicsMaterialSetProperty(materials[2],
                              TRANSPORT_XSECTIONS,
                              SIMPLEXS1,num_groups,0.1*20,0.8)

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
pquad0 = chiCreateProductQuadrature(GAUSS_LEGENDRE_CHEBYSHEV,12, 2)
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
QOI_value_sum = 0.0
for g=0,num_groups-1 do
    ff = chiGetFieldFunctionHandleByName("phi_g"..
            string.format("%03d",g).."_m"..string.format("%02d",0))
    ffi1 = chiFFInterpolationCreate(VOLUME)
    curffi = ffi1
    chiFFInterpolationSetProperty(curffi,OPERATION,OP_SUM)
    chiFFInterpolationSetProperty(curffi,LOGICAL_VOLUME,tvol1)
    chiFFInterpolationSetProperty(curffi,ADD_FIELDFUNCTION,ff)

    chiFFInterpolationInitialize(curffi)
    chiFFInterpolationExecute(curffi)
    QOI_value = chiFFInterpolationGetValue(curffi)

    chiLog(LOG_0,string.format("QOI-value["..tostring(g).."]= %.5e", QOI_value))

    QOI_value_sum = QOI_value_sum + QOI_value
end
chiLog(LOG_0,string.format("QOI-value[sum]= %.5e", QOI_value_sum))

--############################################### Exports
if master_export == nil then
    chiExportMultiFieldFunctionToVTK({ff_m0, ff_m1, ff_m2},"ZPhi_LBS")
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
