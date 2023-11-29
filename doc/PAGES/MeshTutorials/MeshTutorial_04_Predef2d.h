/** \page MeshTutorial_04 Mesh Tutorial 4: Predefined 2D Mesh Setup

Consider an example of using a predeveloped mesh as shown in Figure 1 below.

\image html Meshing/SquareMesh2x2.png "Figure 1 - Mesh spanning from [-1,-1] to [1,1]" width=200px


The first step in this process is to create a chi_mesh::MeshHandler for all
subsequent meshing operations. This will also serve as a global subspace for
other meshing objects that will be added.

\code
MeshHandlerCreate()
\endcode

We next import a surface mesh (chi_mesh::SurfaceMesh) that will define our
 problem. **Solvers will directly operate
on predefined meshes**.

\code
newSurfMesh = SurfaceMeshCreate();
SurfaceMeshImportFromOBJFile(newSurfMesh,"CHI_RESOURCES/TestObjects/SquareMesh2x2.obj")
\endcode

The next step is to break up the boundaries of this mesh into edges that we can
assign to boundaries. We first collect *loops* of edges, which are lists of
 connected edges. We then split these loops by angle into further loops. Finally
 we assign each fine-grained loop to a chi_mesh::LineMesh object that can be
 added to a chi_mesh::Boundary type.

 \code
line_mesh = {};
line_mesh_count = 0;
--
for k=1,loop_count do
  split_loops,split_count = chiEdgeLoopSplitByAngle(loops,k-1);
  for m=1,split_count do
    line_mesh_count = line_mesh_count + 1;
    line_mesh[line_mesh_count] = chiLineMeshCreateFromLoop(split_loops,m-1);
  end
--
end
\endcode

Next we create a chi_mesh::Region object to which we are going to add boundaries.
We can add a boundary for each chi_mesh::LineMesh and each chi_mesh::SurfaceMesh.

\code
region1 = chiRegionCreate()
chiRegionAddSurfaceBoundary(region1,newSurfMesh);
for k=1,line_mesh_count do
  chiRegionAddLineBoundary(region1,line_mesh[k]);
end
\endcode

We now have everything we need to do a Predefined 2D mesh operation. To do this
we create a chi_mesh::SurfaceMesher with the type *SURFACEMESHER_PREDEFINED* as
 well as a chi_mesh::VolumeMesher with the type *VOLUMEMESHER_PREDEFINED2D*.

\code
SurfaceMesherCreate(SURFACEMESHER_PREDEFINED);
VolumeMesherCreate(VOLUMEMESHER_PREDEFINED2D);
--
SurfaceMesherExecute();
VolumeMesherExecute()
\endcode

The *SURFACEMESHER_PREDEFINED* is just a pass-through mesher, in contrast to the
 *SURFACEMESHER_DELAUNAY* type which will re-mesh surfaces. The
 *VOLUMEMESHER_PREDEFINED2D* similarly does not create new mesh elements but at
 least this instantiates the notion of nodes and cells.

### Complete code segment

\code
MeshHandlerCreate()
--
newSurfMesh = SurfaceMeshCreate();
SurfaceMeshImportFromOBJFile(newSurfMesh,"CHI_RESOURCES/TestObjects/SquareMesh2x2.obj")
loops,loop_count = SurfaceMeshGetEdgeLoops(newSurfMesh)
--
--
line_mesh = {};
line_mesh_count = 0;
--
for k=1,loop_count do
  split_loops,split_count = chiEdgeLoopSplitByAngle(loops,k-1);
  for m=1,split_count do
    line_mesh_count = line_mesh_count + 1;
    line_mesh[line_mesh_count] = chiLineMeshCreateFromLoop(split_loops,m-1);
  end
--
end
--
region1 = chiRegionCreate()
chiRegionAddSurfaceBoundary(region1,newSurfMesh);
for k=1,line_mesh_count do
  chiRegionAddLineBoundary(region1,line_mesh[k]);
end
--
SurfaceMesherCreate(SURFACEMESHER_PREDEFINED);
VolumeMesherCreate(VOLUMEMESHER_PREDEFINED2D);
--
SurfaceMesherExecute();
VolumeMesherExecute()
\endcode

*/
