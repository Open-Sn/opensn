/** \page Tutorial03 Tutorial 3: Basic transport simulation

Let us tackle a basic transport simulation on a 3D mesh\n

\f[
 \hat{\Omega} \nabla \Psi(\hat{\Omega}) + \sigma_t \Psi (\hat{\Omega}) =
 \int_{\hat{\Omega}} \sigma_s (\hat{\Omega}' \to \hat{\Omega})
 \Psi (\hat{\Omega}').d\hat{\Omega}' + q(\hat{\Omega})
\f]

## Step 1 - Create your input file, open it in a text editor

 As before, go to the directory of your
 choosing and make an input file of your choosing.

## Step 2 - Create a mesh

As was done in Tutorials 1 and 2 we create the 3D mesh as follows:

\code
MeshHandlerCreate()

nodes={}
N=32
ds=2.0/N
for i=0,N do
    nodes[i+1] = -1.0 + i*ds
end
surf_mesh,region1 = MeshCreateUnpartitioned3DOrthoMesh(nodes,nodes,nodes)

chiVolumeMesherExecute();
\endcode

Next we set the material IDs:

\code
-- Set Material IDs
vol0 = LogicalVolumeCreate(RPP,-1000,1000,-1000,1000,-1000,1000)
VolumeMesherSetProperty(MATID_FROMLOGICAL,vol0,0)
\endcode

## Step 3 - Add a material with a transport cross-section

Similar to the diffusion tutorial we have to create a material, then add a
 property to it, and then to set the property.

\code
material0 = PhysicsAddMaterial("Test Material");

PhysicsMaterialAddProperty(material0,TRANSPORT_XSECTIONS)
num_groups = 1
PhysicsMaterialSetProperty(material0,
                              TRANSPORT_XSECTIONS,
                              SIMPLEXS1,
                              num_groups,     --Num grps
                              1.0,            --Sigma_t
                              0.2)            --Scattering ratio
\endcode

The general material property TRANSPORT_XSECTIONS is used for
 non-fission transport materials. Simple cross-sections can be specified by
 using the operation SIMPLEXS1 which expects the number of groups,
 the total cross-section, and then the scattering ratio.

## Step 4 - Add Transport physics

\code
--############################################### Setup Physics
phys1 = LBSCreateSolver()
SolverAddRegion(phys1,region1)

for k=1,num_groups do
    LBSCreateGroup(phys1)
end

pquad = CreateProductQuadrature(GAUSS_LEGENDRE_CHEBYSHEV,4,4)

--========== Groupset def
gs0 = LBSCreateGroupset(phys1)
LBSGroupsetAddGroups(phys1,gs0,0,num_groups-1)
LBSGroupsetSetQuadrature(phys1,gs0,pquad)
LBSGroupsetSetAngleAggregationType(phys1,gs0,LBSGroupset.ANGLE_AGG_SINGLE)
LBSGroupsetSetIterativeMethod(NPT_GMRES_CYCLES)

--========== Boundary conditions
bsrc = {}
for k=1,num_groups do
    bsrc[k] = 0.0
end
bsrc[1] = 0.5
LBSSetProperty(phys1,BOUNDARY_CONDITION,
                  YMIN,LBSBoundaryTypes.INCIDENT_ISOTROPIC,bsrc);

--========== Solvers
LBSSetProperty(phys1,DISCRETIZATION_METHOD,PWLD1D)
LBSSetProperty(phys1,SCATTERING_ORDER,0)
\endcode

 A transport solver is invoked by using a call to LBSCreateSolver().
 This creates a derived object based on a base physics solver so the
 mesh region gets added to the
 solver using the generic call SolverAddRegion(). Past this point we need
 to create the single required group with LBSCreateGroup(), although we put
 this in a loop for in-case we want to increase the number of groups, and then a
 quadrature rule for integration of the angular fluxes. Since we are dealing
 with a 3D simulation we will be integrating over \f$ \theta \f$, the polar
 angle, and \f$ \varphi \f$, the azimuthal angle. A quadrature with favorable
 parallel properties is the Gauss-Legendre-Chebyshev quadrature. We create this
 quadrature with a call to
 CreateProductQuadrature() and specifying GAUSS_LEGENDRE_CHEBYSHEV as the rule
 and we then specify 2 azimuthal angles per octant and 2 polar angles per octant.

 The next step in the process is to assign a group-set. Group-sets are very
 useful aggregation features in higher dimension simulations but here we
 only have a single groupset. The group-set is created with a call to
 LBSCreateGroupset(). Next we add groups to the group-set using a range,
 however, since we only have one group here the range will be 0 to 0. The
 final piece of a groupset is to add a quadrature which is achieved with a
 call to LBSGroupsetSetQuadrature().


## Step 5 - Initialize and Execute

\code
LBSInitialize(phys1)
LBSExecute(phys1)
\endcode

This should be intuitive.

## Step 6 - Add output

\code
fflist,count = LBSGetScalarFieldFunctionList(phys1)

cline = FFInterpolationCreate(LINE)
FFInterpolationSetProperty(cline,LINE_FIRSTPOINT,0.0,-1.0,-1.0)
FFInterpolationSetProperty(cline,LINE_SECONDPOINT,0.0, 1.0,1.0)
FFInterpolationSetProperty(cline,LINE_NUMBEROFPOINTS, 50)

FFInterpolationSetProperty(cline,ADD_FIELDFUNCTION,fflist[1])


FFInterpolationInitialize(cline)
FFInterpolationExecute(cline)
FFInterpolationExportPython(cline)

if (location_id == 0) then
    local handle = io.popen("python ZLFFI00.py")
end
\endcode

Instead of using a SLICE interpolator we instead opt to use a LINE interpolator.
 A line-interpolator needs two points as the initial and final point (i.e.
 LINE_FIRSTPOINT and LINE_SECONDPOINT) as well as the total number of
 points LINE_NUMBEROFPOINTS. Finally add the first available scalar field
 field function to the interpolator.

 The output should look as follows:

  \image html "Physics/Tut3Output.png"  "Figure 1 - Output of the 1D simulation" width=350px

## Fully commented code

\code
--############################################### Setup mesh
MeshHandlerCreate()

nodes={}
N=32
ds=2.0/N
for i=0,N do
    nodes[i+1] = -1.0 + i*ds
end
surf_mesh,region1 = MeshCreateUnpartitioned3DOrthoMesh(nodes,nodes,nodes)

--chiSurfaceMesherSetProperty(PARTITION_X,2)
--chiSurfaceMesherSetProperty(PARTITION_Y,2)
--chiSurfaceMesherSetProperty(CUT_X,0.0)
--chiSurfaceMesherSetProperty(CUT_Y,0.0)

chiVolumeMesherExecute();

-- Set Material IDs
vol0 = LogicalVolumeCreate(RPP,-1000,1000,-1000,1000,-1000,1000)
VolumeMesherSetProperty(MATID_FROMLOGICAL,vol0,0)

--############################################### Add material
material0 = PhysicsAddMaterial("Test Material");

PhysicsMaterialAddProperty(material0,TRANSPORT_XSECTIONS)
num_groups = 1
PhysicsMaterialSetProperty(material0,
                              TRANSPORT_XSECTIONS,
                              SIMPLEXS1,
                              num_groups,     --Num grps
                              1.0,   --Sigma_t
                              0.2)   --Scattering ratio

--############################################### Setup Physics
phys1 = LBSCreateSolver()
SolverAddRegion(phys1,region1)

for k=1,num_groups do
    LBSCreateGroup(phys1)
end

pquad = CreateProductQuadrature(GAUSS_LEGENDRE_CHEBYSHEV,2,2)

--========== Groupset def
gs0 = LBSCreateGroupset(phys1)
LBSGroupsetAddGroups(phys1,gs0,0,num_groups-1)
LBSGroupsetSetQuadrature(phys1,gs0,pquad)
LBSGroupsetSetAngleAggregationType(phys1,gs0,LBSGroupset.ANGLE_AGG_SINGLE)
LBSGroupsetSetIterativeMethod(phys1,gs0,NPT_GMRES_CYCLES)

--========== Boundary conditions
bsrc = {}
for k=1,num_groups do
    bsrc[k] = 0.0
end
bsrc[1] = 0.5
LBSSetProperty(phys1,BOUNDARY_CONDITION,
                  YMIN,LBSBoundaryTypes.INCIDENT_ISOTROPIC,bsrc);

--========== Solvers
LBSSetProperty(phys1,DISCRETIZATION_METHOD,PWLD)
LBSSetProperty(phys1,SCATTERING_ORDER,0)

LBSInitialize(phys1)
LBSExecute(phys1)

--############################################### Setup Output
fflist,count = LBSGetScalarFieldFunctionList(phys1)

cline = FFInterpolationCreate(LINE)
FFInterpolationSetProperty(cline,LINE_FIRSTPOINT,0.0,-1.0,-1.0)
FFInterpolationSetProperty(cline,LINE_SECONDPOINT,0.0, 1.0,1.0)
FFInterpolationSetProperty(cline,LINE_NUMBEROFPOINTS, 50)

FFInterpolationSetProperty(cline,ADD_FIELDFUNCTION,fflist[1])


FFInterpolationInitialize(cline)
FFInterpolationExecute(cline)
FFInterpolationExportPython(cline)

ExportFieldFunctionToVTK(fflist[1],"Tutorial3Output","Phi")

if (location_id == 0) then
    local handle = io.popen("python ZLFFI00.py")
end
\endcode


*/
