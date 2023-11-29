#pragma once

#include "framework/mesh/unpartitioned_mesh/unpartitioned_mesh.h"

#include "framework/lua.h"

namespace opensnlua
{

/**Creates an empty unpartitioned mesh. An empty unpartitioned mesh
 * is meant to be manipulated with calls to UnpartitionedMeshUploadVertex()
 * and UnpartitionedMeshUploadCell(). It essentially supports building a mesh
 * manually.
 *
 * ##_
 *
 * ###Example
 * Example usage
 * \code
 * umesh = CreateEmptyUnpartitionedMesh()
 * \endcode
 *
 * \ingroup LuaUnpartitionedMesh
 */
int CreateEmptyUnpartitionedMesh(lua_State* L);

/**Destroy an unpartitioned mesh. This routine should be called for
 * memory sensitive simulations because each process will have a full
 * copy of this data.
 *
 * \param handle int Handle to mesh.
 *
 * ##_
 *
 * ###Example
 * Example usage
 * \code
 * DestroyUnpartitionedMesh(umesh)
 * \endcode
 *
 * \ingroup LuaUnpartitionedMesh
 */
int DestroyUnpartitionedMesh(lua_State* L);

/**Creates an unpartitioned mesh from VTK Unstructured mesh files.
 *
 * \param file_name char Filename of the .vtu file.
 * \param field char Name of the cell data field from which to read
 *                   material and boundary identifiers (optional).
 *
 * \ingroup LuaUnpartitionedMesh
 *
 * ##_
 *
 * ### Example
 * An example mesh creation below:
 * \code
 * MeshHandlerCreate()
 *
 * umesh = UnpartitionedMeshFromVTU("ZMeshTest_0.vtu")
 *
 * SurfaceMesherCreate(SURFACEMESHER_PREDEFINED)
 * chiVolumeMesherCreate(VOLUMEMESHER_UNPARTITIONED, umesh)
 *
 * SurfaceMesherExecute()
 * chiVolumeMesherExecute()
 * \endcode
 *
 *
 * \return Handle A handle to the newly created UnpartitionedMesh
 */
int UnpartitionedMeshFromVTU(lua_State* L);

/**Creates an unpartitioned mesh from VTK Partitioned Unstructured mesh files
 * (.pvtu).
 *
 * \param file_name char Filename of the .vtu file.
 * \param field char Name of the cell data field from which to read
 *                   material and boundary identifiers (optional).
 *
 * \ingroup LuaUnpartitionedMesh
 *
 * ##_
 *
 * ### Example
 * An example mesh creation below:
 * \code
 * MeshHandlerCreate()
 *
 * umesh = UnpartitionedMeshFromPVTU("ZMeshTest_0.vtu")
 *
 * SurfaceMesherCreate(SURFACEMESHER_PREDEFINED)
 * chiVolumeMesherCreate(VOLUMEMESHER_UNPARTITIONED, umesh)
 *
 * SurfaceMesherExecute()
 * chiVolumeMesherExecute()
 * \endcode
 *
 *
 * \return Handle A handle to the newly created UnpartitionedMesh
 */
int UnpartitionedMeshFromPVTU(lua_State* L);

/**Creates an unpartitioned mesh from starccm+ exported
 * Ensight Gold mesh files.
 *
 * \param file_name char Filename of the .case file.
 * \param scale float Scale to apply to the mesh
 *
 * \ingroup LuaUnpartitionedMesh
 *
 * ##_
 *
 * ### Example
 * An example mesh creation below:
 * \code
 * MeshHandlerCreate()
 *
 * umesh = UnpartitionedMeshFromEnsightGold("resources/TestObjects/Sphere.case")
 *
 * SurfaceMesherCreate(SURFACEMESHER_PREDEFINED)
 * chiVolumeMesherCreate(VOLUMEMESHER_UNPARTITIONED, umesh)
 *
 * SurfaceMesherExecute()
 * chiVolumeMesherExecute()
 * \endcode
 *
 * \return Handle A handle to the newly created UnpartitionedMesh
 */
int UnpartitionedMeshFromEnsightGold(lua_State* L);

/**Creates an unpartitioned mesh from a wavefront .obj file.
 *
 * \param file_name char Filename of the .case file.
 *
 * \ingroup LuaUnpartitionedMesh
 *
 * ##_
 *
 * ### Example
 * An example mesh creation below:
 * \code
 * MeshHandlerCreate()
 *
 * umesh =
 * UnpartitionedMeshFromWavefrontOBJ("resources/TestObjects/TriangleMesh2x2.obj")
 *
 * SurfaceMesherCreate(SURFACEMESHER_PREDEFINED)
 * chiVolumeMesherCreate(VOLUMEMESHER_UNPARTITIONED, umesh)
 *
 * SurfaceMesherExecute()
 * chiVolumeMesherExecute()
 * \endcode
 *
 * \return Handle A handle to the newly created UnpartitionedMesh
 */
int UnpartitionedMeshFromWavefrontOBJ(lua_State* L);

/**Creates an unpartitioned mesh from a .msh file.
 *
 * \param file_name char Filename of the .msh file.
 *
 * \ingroup LuaUnpartitionedMesh
 *
 * ##_
 *
 * ### Example
 * An example mesh creation below:
 * \code
 * MeshHandlerCreate()
 *
 * umesh = UnpartitionedMeshFromMshFormat("File.msh")
 *
 * SurfaceMesherCreate(SURFACEMESHER_PREDEFINED)
 * chiVolumeMesherCreate(VOLUMEMESHER_UNPARTITIONED, umesh)
 *
 * SurfaceMesherExecute()
 * chiVolumeMesherExecute()
 * \endcode
 *
 *
 * \return Handle A handle to the newly created UnpartitionedMesh
 */
int UnpartitionedMeshFromMshFormat(lua_State* L);

/**Creates an unpartitioned mesh from ExodusII format.
 *
 * \param file_name char Filename of the .case file.
 * \param scale float Scale to apply to the mesh
 *
 * \ingroup LuaUnpartitionedMesh
 *
 * ##_
 *
 * ### Example
 * An example mesh creation below:
 * \code
 * MeshHandlerCreate()
 *
 * umesh = UnpartitionedMeshFromExodusII("resources/TestObjects/Mesh.e")
 *
 * SurfaceMesherCreate(SURFACEMESHER_PREDEFINED)
 * chiVolumeMesherCreate(VOLUMEMESHER_UNPARTITIONED, umesh)
 *
 * SurfaceMesherExecute()
 * chiVolumeMesherExecute()
 * \endcode
 *
 * \return Handle A handle to the newly created UnpartitionedMesh
 */
int UnpartitionedMeshFromExodusII(lua_State* L);

/**Uploads a vertex.
 *
 * \param handle int Handle to mesh.
 * \param x      double x-coordinate.
 * \param y      double y-coordinate.
 * \param z      double z-coordinate.
 *
 * ## _
 *
 * ###Example
 * Example usage
 * \code
 * UnpartitionedMeshUploadVertex(umesh, 0, 0, 0)
 * UnpartitionedMeshUploadVertex(umesh, 1, 0, 0)
 * UnpartitionedMeshUploadVertex(umesh, 1, 1, 0)
 * UnpartitionedMeshUploadVertex(umesh, 0, 1, 0)
 * \endcode
 *
 * \ingroup LuaUnpartitionedMesh
 */
int UnpartitionedMeshUploadVertex(lua_State* L);

/**Uploads a cell
 *
 * \param handle int Handle to mesh.
 * \param cell_table lua_table A Lua-table containing fields of data. See
 *  cell_table below.
 *
 * ## _
 *
 * ###cell_table
 * A lua-table with the following fields:
 * - `type` <I>string</I>. The value of this field contains the cell's primary
 * type. The value can be "SLAB", "POLYGON" or "POLYHEDRON".
 * - `sub_type` <I>string</I>. The value of this field constains the cell's
 *  secondary type. The value can be "SLAB", "POLYGON", "TRIANGLE",
 *  "QUADRILATERAL", "POLYHEDRON", "TETRAHEDRON" or "HEXAHEDRON".
 * - `num_faces` <I>int</I>. The value of this field represent the number of faces
 *  specified for this cell. Each face is contained in a field "faceX" where X is
 *  the face index.
 * - `material_id` <I>int</I>. (Optional) The value of this field holds a material
 *  identifier. If not provided, will be defaulted to -1.
 * - `faceX` <I>table</I>. A field holding a lua-table containing the vertex-ids
 *  of face X. There must be `num_faces` of these fields, i.e., "face0", "face1",
 * etc.
 *
 * ###Example
 * Example usage
 * \code
 * UnpartitionedMeshUploadVertex(umesh, 0, 0, 0)
 * UnpartitionedMeshUploadVertex(umesh, 1, 0, 0)
 * UnpartitionedMeshUploadVertex(umesh, 1, 1, 0)
 * UnpartitionedMeshUploadVertex(umesh, 0, 1, 0)
 *
 * UnpartitionedMeshUploadVertex(umesh, 0, 0, 1)
 * UnpartitionedMeshUploadVertex(umesh, 1, 0, 1)
 * UnpartitionedMeshUploadVertex(umesh, 1, 1, 1)
 * UnpartitionedMeshUploadVertex(umesh, 0, 1, 1)
 *
 * cell = {}
 * cell.type        = "POLYHEDRON"
 * cell.sub_type    = "HEXAHEDRON"
 * cell.num_faces   = 6
 * cell.material_id = 0
 * cell.face0 = {1,2,6,5}
 * cell.face1 = {0,4,7,3}
 * cell.face2 = {2,3,7,6}
 * cell.face3 = {0,1,5,4}
 * cell.face4 = {4,5,6,7}
 * cell.face5 = {0,3,2,1}
 *
 * UnpartitionedMeshUploadCell(umesh, cell, true)
 * UnpartitionedMeshFinalizeEmpty(umesh)
 * \endcode
 *
 * \ingroup LuaUnpartitionedMesh
 */
int UnpartitionedMeshUploadCell(lua_State* L);

/**Finalizes a mesh. This usually involves computing centroids and
 * establishing connectivity.
 *
 * \param handle int Handle to mesh.
 *
 * ## _
 *
 * ###Example
 * Example usage
 * \code
 * UnpartitionedMeshUploadVertex(umesh, 0, 0, 0)
 * UnpartitionedMeshUploadVertex(umesh, 1, 0, 0)
 * UnpartitionedMeshUploadVertex(umesh, 1, 1, 0)
 * UnpartitionedMeshUploadVertex(umesh, 0, 1, 0)
 *
 * UnpartitionedMeshUploadVertex(umesh, 0, 0, 1)
 * UnpartitionedMeshUploadVertex(umesh, 1, 0, 1)
 * UnpartitionedMeshUploadVertex(umesh, 1, 1, 1)
 * UnpartitionedMeshUploadVertex(umesh, 0, 1, 1)
 *
 * cell = {}
 * cell.type        = "POLYHEDRON"
 * cell.sub_type    = "HEXAHEDRON"
 * cell.num_faces   = 6
 * cell.material_id = 0
 * cell.face0 = {1,2,6,5}
 * cell.face1 = {0,4,7,3}
 * cell.face2 = {2,3,7,6}
 * cell.face3 = {0,1,5,4}
 * cell.face4 = {4,5,6,7}
 * cell.face5 = {0,3,2,1}
 *
 * UnpartitionedMeshUploadCell(umesh, cell, true)
 * UnpartitionedMeshFinalizeEmpty(umesh)
 * \endcode
 *
 * \ingroup LuaUnpartitionedMesh
 */
int UnpartitionedMeshFinalizeEmpty(lua_State* L);

} // namespace opensnlua
