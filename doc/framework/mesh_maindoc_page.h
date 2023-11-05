/** \defgroup LuaMesh A Meshing


## Mesh Handling

 Meshing in ChiTech is made available via the chi_mesh namespace. All instances
 of entities are stored within a handler (chi_mesh::MeshHandler). A mesh handler
 is created with a call to chiMeshHandlerCreate().

 \code
 chiMeshHandlerCreate()
 \endcode

## Predefined 2D Mesh Setup Example

\code
chiMeshHandlerCreate()
--
newSurfMesh = chiSurfaceMeshCreate();
chiSurfaceMeshImportFromOBJFile(newSurfMesh,"CHI_RESOURCES/TestObjects/SquareMesh2x2.obj")
loops,loop_count = chiSurfaceMeshGetEdgeLoops(newSurfMesh)
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
\endcode

*/
