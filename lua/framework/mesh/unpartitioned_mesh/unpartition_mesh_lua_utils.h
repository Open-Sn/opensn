#pragma once

#include "framework/mesh/unpartitioned_mesh/unpartitioned_mesh.h"

#include "framework/lua.h"

namespace opensnlua
{

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
 * VolumeMesherCreate(VOLUMEMESHER_UNPARTITIONED, umesh)
 *
 * SurfaceMesherExecute()
 * VolumeMesherExecute()
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
 * VolumeMesherCreate(VOLUMEMESHER_UNPARTITIONED, umesh)
 *
 * SurfaceMesherExecute()
 * VolumeMesherExecute()
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
 * VolumeMesherCreate(VOLUMEMESHER_UNPARTITIONED, umesh)
 *
 * SurfaceMesherExecute()
 * VolumeMesherExecute()
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
 * VolumeMesherCreate(VOLUMEMESHER_UNPARTITIONED, umesh)
 *
 * SurfaceMesherExecute()
 * VolumeMesherExecute()
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
 * VolumeMesherCreate(VOLUMEMESHER_UNPARTITIONED, umesh)
 *
 * SurfaceMesherExecute()
 * VolumeMesherExecute()
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
 * VolumeMesherCreate(VOLUMEMESHER_UNPARTITIONED, umesh)
 *
 * SurfaceMesherExecute()
 * VolumeMesherExecute()
 * \endcode
 *
 * \return Handle A handle to the newly created UnpartitionedMesh
 */
int UnpartitionedMeshFromExodusII(lua_State* L);

} // namespace opensnlua
