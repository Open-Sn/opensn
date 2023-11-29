/** \page Tutorial01 Tutorial 1: Basic Diffusion Simulation

 Let us tackle the basic diffusion equation as applied to steady state heat
 transfer on a 3D orthogonal mesh.\n
 \f[
 -\nabla \cdot k \nabla T = q
 \f]

 The purpose of this tutorial is to understand the basics of the
 input paradigm which is completely different from the traditional
 input-card based method.
 ChiTech uses an imbedded Lua console. For a history of Lua please see its
 wikipedia page
 <a href="https://en.wikipedia.org/wiki/Lua_(programming_language)">here</a>.
 This paradigm allows the user to control many aspects of the program flow
 as well to be more versatile when creating inputs, i.e. reading files, for
  loops to create multiple entities etc. So let us get started.

 ## Step 1 - Create your input file, open it in a text editor

 Go to the directory of your
 choosing and make an input file of your choosing. The input file does not
 require a specific .lua or .inp extension although if you use .lua then
 there are many
 text editors that will be able to provide syntax highlighting.

## Step 2 - Create a mesh handler

\code
chiMeshHandlerCreate()
\endcode

This function call creates a new mesh handler and pushes it to a global
stack. It also make the newly created handler the "current" meaning all
 mesh operations will operate on this handler.
 The concept of creating, reading or manipulating computational meshes is
centered on operating on a specific chi_mesh::MeshHandler object. Physics/math
 objects can then freely operate on meshes by simply specifying which handler
 to use.

## Step 3 - Create a 3D mesh

 For this we will create a simple orthogonal grid of 32x32x32 cells
 (32,768 total) as shown in Figure 1.


 \image html "Meshing/CubeMeshTut1.png"  "Figure 1 - Mesh spanning from [-1,-1,-1] to [1,1,1]" width=350px

The first step of this process is to define an array of nodes starting from -1.0
with a distance \f$ ds \f$ between each of them. In lua we can create this array
programatically:
\code
nodes={}
N=32
ds=2.0/N
for i=0,N do
    nodes[i+1] = i*ds
end
\endcode

Note here that lua arrays cannot use `nodes[0]` since arrays in lua are not
zero based. The next step is call chiMeshCreateUnpartitioned3DOrthoMesh with
the node arrangement along the X-direction, Y-direction and Z-direction,
which in our case is the same,

\code
surf_mesh, region1 = chiMeshCreateUnpartitioned3DOrthoMesh(nodes,nodes,nodes)
\endcode

Some work has been encapsulated behind the scenes here. Firstly, a 2D surface
mesh was created and assigned to a region.
Secondly a surface-mesher (chi_mesh::SurfaceMesher)
 and a volume mesher (chi_mesh::VolumeMesher) has been
assigned to the mesh handler. The volume mesher is not executed because the
macro allows the user to first assign parallel partitioning.

## Step 4 -  Executing the mesher
The underlying mesher for a 3D mesh is an extruded mesh from 2D. This mesh
 has some partitioning properties that can be set to allow parallel simulations,
 however, we will just focus on serial simulations for now.
 The next step to complete meshing is to execute the volume mesher. Since we
 are not going to set any partitioning parameters we will just execute it now.

\code
chiVolumeMesherExecute();
\endcode

## Step 5 - Assign material- and boundary IDs and

Materials ID's can conveniently be specified using logical volumes. There are
many options for logical volumes ranging from primitive parametric surfaces
 to non-convex user generated surfaces using surface meshes
 (see chiLogicalVolumeCreate). For this tutorial
 we will use a rectangular paralellipiped (RPP or brick) as follows.

\code
material = chiPhysicsAddMaterial("Test Material");

vol0 = chiLogicalVolumeCreate(RPP,-1000,1000,-1000,1000,-1000,1000)
VolumeMesherSetProperty(MATID_FROMLOGICAL,vol0,material)
\endcode

We first create a material using the chiPhysicsAddMaterial() function.
 The handle of this material will essentially be zero but it is not technically
 required by the VolumeMesherSetProperty() function since this function
 operates on the integer supplied.
This material is added to the physics environment and therefore has scope over
all mesh handlers and all physics entities. We then create a logical volume
 using the function chiLogicalVolumeCreate() with arguments RPP,
 specifying that we will be using a Rectangular Parallelipiped and then a
 series of dimensions specifying xmin-xmax-ymin-ymax-zmin-zmax. Every cell
 centroid within these dimensions will be flagged as being "within" the logical
 volume.

Actually setting the mesh's material id's is a utility facilitated by a
 volume mesher property and hence we call VolumeMesherSetProperty()
 with a property index MATID_FROMLOGICAL. Next we provided a handle to the
 logical volume (vol0) and the desired material id. Logical volumes are very
 diverse and their uses are discussed elsewhere. There is an additional utility,
 chiVolumeMesherSetMatIDToAll(), which will set all cell-material-ids without
 requiring a logical volume, however, in this tutorial we opted to show the
 logical volume route.

### Boundary IDs
By default all boundaries are unassigned (i.e., -1). There are two utilities
 that set the boundary id's the first of which is essentially identical to how
 material ids are set, VolumeMesherSetProperty(), but this time the property
 index is BNDRYID_FROMLOGICAL. The second way is to use
 chiVolumeMesherSetupOrthogonalBoundaries() which requires no arguments and will
 asssign standard indices to a boundary if it is aligned with the orthogonal
 cartesian directions.

The culmination of this step is all done within a physics agnostic framework.
The user can even export the mesh for visualization using the function
 chiMeshHandlerExportMeshToObj().

### Boundary IDs


## Step 5 - Adding material properties

Now that the cells have been assigned a material id we need to add
 properties to the material conducive to a diffusion simulation. For a heat
 transfer diffusion simulation we will need to know the thermal conductivity
 "k" and the volumetric source strength "q". Both of these can be simple
 scalar values.

\code
chiPhysicsMaterialAddProperty(material,SCALAR_VALUE,"k")
chiPhysicsMaterialSetProperty(material,"k",SINGLE_VALUE,1.0)

chiPhysicsMaterialAddProperty(material,SCALAR_VALUE,"q")
chiPhysicsMaterialSetProperty(material,"q",SINGLE_VALUE,1.0)
\endcode

In this code we created the material properties using the function
 chiPhysicsMaterialAddProperty() which requires a handle to the reference
 material, the property type (SCALAR_VALUE), and a name for the property.
 Unlike extrusion layers the property name can be used in further calls to refer
 to the specific property.

Material property values are set using the function
 chiPhysicsMaterialSetProperty() which again expects a handle to the
 reference material, then either a material property id or name (in this case
 name), then an operation index and value(s). For this case we used an operation
 index SINGLE_VALUE which is the only operation supported by SCALAR_VALUE. In
 future the user can specify, as an example, temperature dependent values which
 will support the operation FROM_TABLE, but that is a topic for a different time.

## Step 6 - Setup the diffusion physics

The following sequence of function calls completely define the diffusion solver.

\code
phys1 = chiDiffusionCreateSolver()
chiSolverSetBasicOption(phys1,"discretization_method","PWLC")
chiSolverSetBasicOption(phys1,"residual_tolerance",1.0e-6)
\endcode

We first create the diffusion solver with a call to
 chiDiffusionCreateSolver(). This creates the solver and pushes it onto
 the physics handler. The function returns the handle.

Next we can set numerous diffusion solver properties which can comprehensively
 be viewed in its specific documentation (chiDiffusionSetProperty()).

## Step 7 - Initialize and Solve

The final step of this process is to initialize and execute the diffusion solver.

\code
chiDiffusionInitialize(phys1)
chiDiffusionExecute(phys1)
\endcode

## Step 8 - Post-processing
The execution of a Chi-Tech physics module culminates in the creation of one or
 more field functions. These field functions should all be populated for use
 during the initialization phase of the solver. Users can get handles to all
 the field functions by using the chiGetFieldFunctionList call. Thereafter, the
 field function can be exported to VTK-format which can be read by Paraview.

\code
fflist,count = chiGetFieldFunctionList(phys1)
chiExportFieldFunctionToVTK(fflist[1],"Tutorial1Output","Temperature")
\endcode



## Complete input file

Here is the complete input file with comments

\code
--############################################### Setup mesh
chiMeshHandlerCreate()

nodes={}
N=32
ds=2.0/N
for i=0,N do
    nodes[i+1] = -1.0 + i*ds
end
surf_mesh,region1 = chiMeshCreateUnpartitioned3DOrthoMesh(nodes,nodes,nodes)

chiVolumeMesherExecute();

material = chiPhysicsAddMaterial("Test Material");

-- Set Material IDs
vol0 = chiLogicalVolumeCreate(RPP,-1000,1000,-1000,1000,-1000,1000)
VolumeMesherSetProperty(MATID_FROMLOGICAL,vol0,material)

MeshHandlerExportMeshToVTK("Mesh")
--############################################### Add material properties


-- Set material properties
chiPhysicsMaterialAddProperty(material,SCALAR_VALUE,"k")
chiPhysicsMaterialSetProperty(material,"k",SINGLE_VALUE,1.0)

chiPhysicsMaterialAddProperty(material,SCALAR_VALUE,"q")
chiPhysicsMaterialSetProperty(material,"q",SINGLE_VALUE,1.0)


--############################################### Setup Physics
phys1 = chiDiffusionCreateSolver()
chiSolverSetBasicOption(phys1,"discretization_method","PWLC");
chiSolverSetBasicOption(phys1,"residual_tolerance",1.0e-6)

--############################################### Initialize and
--                                                Execute Solver
chiDiffusionInitialize(phys1)
chiDiffusionExecute(phys1)

----############################################### Visualize the field function
fflist,count = chiGetFieldFunctionList(phys1)
chiExportFieldFunctionToVTK(fflist[1],"Tutorial1Output","Temperature")
\endcode


## Step 9 - Execute the code
Assuming you added the executable to your PATH environment variable, the code
can be executed by typing the executable name followed by the input file path
 (relative or absolute).

\verbatim
ChiTech Tutorial01.lua
\endverbatim

The output produced will look as follows:

\verbatim
[0]  Parsing argument 1 Tutorials/Tutorial01.lua
[0]  2022-05-09 14:01:21 Running ChiTech in batch-mode with 1 processes.
[0]  ChiTech number of arguments supplied: 1
[0]  Computing cell-centroids.
[0]  Done computing cell-centroids.
[0]  Checking cell-center-to-face orientations
[0]  Done checking cell-center-to-face orientations
[0]  00:00:00 Establishing cell connectivity.
[0]  00:00:00 Vertex cell subscriptions complete.
[0]  00:00:00 Surpassing cell 3277 of 32768 (10%)
[0]  00:00:00 Surpassing cell 6554 of 32768 (20%)
[0]  00:00:00 Surpassing cell 9831 of 32768 (30%)
[0]  00:00:00 Surpassing cell 13108 of 32768 (40%)
[0]  00:00:00 Surpassing cell 16384 of 32768 (50%)
[0]  00:00:00 Surpassing cell 19661 of 32768 (60%)
[0]  00:00:00 Surpassing cell 22938 of 32768 (70%)
[0]  00:00:00 Surpassing cell 26215 of 32768 (80%)
[0]  00:00:00 Surpassing cell 29492 of 32768 (90%)
[0]  00:00:00 Surpassing cell 32768 of 32768 (100%)
[0]  00:00:00 Establishing cell boundary connectivity.
[0]  00:00:00 Done establishing cell connectivity.
[0]  00:00:00 VolumeMesherPredefinedUnpartitioned executing. Memory in use = 61 MB
[0]  Computed centroids
[0]  Partitioning mesh with ParMETIS.
[0]  Done partitioning mesh.
[0]  Cells loaded.
[0]  VolumeMesherPredefinedUnpartitioned: Cells created = 32768
[0]  00:00:00 chiVolumeMesherExecute: Volume meshing completed. Memory used = 43.5 MB
[0]  Total process memory used after meshing 104 MB
[0]  00:00:00 Setting material id from logical volume.
[0]  00:00:00 Done setting material id from logical volume. Number of cells modified = 32768.
[0]  Exporting mesh to VTK. 32768
[0]  Done exporting mesh to VTK.
[0]  Solver:DiffusionSolver option:discretization_method set to PWLC.
[0]  Solver:DiffusionSolver option:residual_tolerance set to 1e-06.
[0]
[0]  00:00:00 DiffusionSolver: Initializing Diffusion solver
[0]  Identifying unique boundary-ids.
[0]  Max boundary id identified: 0
[0]  00:00:00 Creating Piecewise Linear Continuous Finite Element spatial discretizaiton.
[0]  00:00:08 Developing nodal ordering.
[0]  00:00:08 Done creating Piecewise Linear Continuous Finite Element spatial discretizaiton.
[0]  DiffusionSolver: Global number of DOFs=35937
[0]  Building sparsity pattern.
[0]  00:00:09 DiffusionSolver: Diffusion Solver initialization time 8.76791
[0]  Setting matrix preallocation.
[0]  00:00:09 DiffusionSolver: Assembling A locally
[0]  00:00:09 DiffusionSolver: Done Assembling A locally
[0]  00:00:09 DiffusionSolver: Communicating matrix assembly
[0]  00:00:09 DiffusionSolver: Assembling A globally
[0]  00:00:09 DiffusionSolver: Diagonal check
[0]  Number of mallocs used = 0
[0]  Number of non-zeros allocated = 912673
[0]  Number of non-zeros used = 759717
[0]  Number of unneeded non-zeros = 152956
[0]  00:00:09 DiffusionSolver: Assembling x and b
[0]  00:00:09 DiffusionSolver: Solving system
[0]  Diffusion iteration    0 - Residual 2.6521952e+02
[0]  Iteration 0 Residual 265.22
[0]  Diffusion iteration    1 - Residual 3.1006314e+01
[0]  Iteration 1 Residual 31.0063
[0]  Diffusion iteration    2 - Residual 4.0467718e+00
[0]  Iteration 2 Residual 4.04677
[0]  Diffusion iteration    3 - Residual 4.9168004e-01
[0]  Iteration 3 Residual 0.49168
[0]  Diffusion iteration    4 - Residual 5.5399485e-02
[0]  Iteration 4 Residual 0.0553995
[0]  Diffusion iteration    5 - Residual 7.3128018e-03
[0]  Iteration 5 Residual 0.0073128
[0]  Diffusion iteration    6 - Residual 1.3610040e-03
[0]  Iteration 6 Residual 0.001361
[0]  Diffusion iteration    7 - Residual 1.7105415e-04
[0]  Iteration 7 Residual 0.000171054
[0]  Diffusion iteration    8 - Residual 2.2641339e-05
[0]  Iteration 8 Residual 2.26413e-05
[0]  Diffusion iteration    9 - Residual 3.3728166e-06
[0]  Iteration 9 Residual 3.37282e-06
[0]  Diffusion iteration   10 - Residual 3.1009104e-07
[0]  Iteration 10 Residual 3.10091e-07
[0]  Convergence reason: KSP_CONVERGED_RTOL
[0]  00:00:09 DiffusionSolver[g=0-0]: Number of iterations =10
[0]  Timing:
[0]  Assembling the matrix: 0.085187
[0]  Solving the system   : 0.075084
[0]  Diffusion Solver execution completed!
[0]  Exporting field function phi to files with base name Tutorial1Output
[0]  Final program time 00:00:10
[0]  2022-05-09 14:01:32 ChiTech finished execution of Tutorials/Tutorial01.lua
\endverbatim

\image html "Meshing/CubeSolutionTut1.png"  "Figure 2 - Solution viewed in Paraview" width=700px


 */
