/** \page Tutorial02 Tutorial 2: Other forms of output
 *
 * Solver modules in ChiTech connect their solution information to a concept
 * called a <B>Field Function</B>. A field function is considered fully defined
 * when it is connected to both a grid (mesh) and a spatial discretization.

 \image html "Physics/FieldFunctionHierarchy.png" "Figure 1 - Hierarchy of field functions"
width=600px

## Step 1 - Make a copy of Tutorial01 input
In the same folder (or any of your choice) make a copy of the input you used for Tutorial 1.
We will be adding some items to this input file.

## Step 2 - Obtain a list of field functions associated with the solver

 \code
 fflist,count = GetFieldFunctionList(phys1)
 \endcode

 The function call GetFieldFunctionList() provides us with two items. A
 lua-table and a count of how many items there are in the table. The items
 in "fflist" are the text names of the field functions. Each solver has its
 own defaults.

## Step 3 - Create a slice interpolator

\code
slice1 = FFInterpolationCreate(SLICE)
FFInterpolationSetProperty(slice1,SLICE_POINT,0.0,0.0,0.0)
FFInterpolationSetProperty(slice1,ADD_FIELDFUNCTION,fflist[1])
\endcode

The first call creates a interpolator of type SLICE. At the time of writing this
 tutorial we support LINE, SLICE and VOLUME. The default orientation of a slice
 interpolator is with the cutting plane's normal pointing in the direction
 of \f$ \hat{k} \f$ and the reference point at (0,0,0). A change in reference
 point is achieved with a call to FFInterpolationSetProperty() with
 a property index SLICE_POINT. The last line here is to add a field function to
 this interpolator. We use the same function but this time with a property index
 ADD_FIELDFUNCTION.

## Step 4 - Initialize and execute
For very complex meshes it might be prudent to perform initialization before
 actually solving the systems, since
 this established the necessary interpolation parameters and allows one to
 execute the interpolator multiple times after that with minimal cost.

\code
FFInterpolationInitialize(slice1)
FFInterpolationExecute(slice1)
FFInterpolationExportPython(slice1)
\endcode

The FFInterpolationInitialize() and FFInterpolationExecute() should be
 intuitive to understand. The last function call here is FFInterpolationExportPython()
 which is a utility to export a slice to a python file. Inside this python
 file there are some default visualization commands but the major utility here
 is that the field function is now represented as python variables so that
 the user can define custom visualizations.

## Step 5 - Run the python file
Being a useful scripting system, the lua console can itself invoke processes.
 Here the default export name of the python file will be "ZPFFI00.py". Its
 actually "ZPFFI0", for Plane-Field-Function-Interpolator-0 but the last digit
 denotes the processor identification in
 such a way that the user merely needs to execute the 0-index and the rest will
 be loaded.

\code
local handle = io.popen("python ZPFFI00.py")
\endcode

The output produced is shown below:

 \image html "Physics/FFOutput.png" "Figure 2 - Output produced by python script" width=600px

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

-- VolumeMesherSetProperty(PARTITION_TYPE,KBA_STYLE_XYZ)
-- VolumeMesherSetKBAPartitioningPxPyPz(2,2,1)
-- VolumeMesherSetKBACutsX({0.0})
-- VolumeMesherSetKBACutsY({0.0})

VolumeMesherExecute();

material = PhysicsAddMaterial("Test Material");

-- Set Material IDs
vol0 = LogicalVolumeCreate(RPP,-1000,1000,-1000,1000,-1000,1000)
mesh.SetMaterialIDFromLogicalVolume(vol0,material)

chiRegionExportMeshToVTK(region1,"Mesh")
--############################################### Add material properties


-- Set material properties
PhysicsMaterialAddProperty(material,SCALAR_VALUE,"k")
PhysicsMaterialSetProperty(material,"k",SINGLE_VALUE,1.0)

PhysicsMaterialAddProperty(material,SCALAR_VALUE,"q")
PhysicsMaterialSetProperty(material,"q",SINGLE_VALUE,1.0)


--############################################### Setup Physics
phys1 = DiffusionCreateSolver();
SolverAddRegion(phys1,region1)
SolverSetBasicOption(phys1,"discretization_method","PWLC");
SolverSetBasicOption(phys1,"residual_tolerance",1.0e-6)

--############################################### Initialize and
--                                                Execute Solver
DiffusionInitialize(phys1)
DiffusionExecute(phys1)

----############################################### Visualize the field function
fflist,count = GetFieldFunctionList(phys1)
ExportFieldFunctionToVTK(fflist[1],"Tutorial1Output","Temperature")

slice1 = FFInterpolationCreate(SLICE)
FFInterpolationSetProperty(slice1,SLICE_POINT,0.0,0.0,0.0)
FFInterpolationSetProperty(slice1,ADD_FIELDFUNCTION,fflist[1])

FFInterpolationInitialize(slice1)
FFInterpolationExecute(slice1)
FFInterpolationExportPython(slice1)

local handle = io.popen("python ZPFFI00.py")
\endcode

 */
