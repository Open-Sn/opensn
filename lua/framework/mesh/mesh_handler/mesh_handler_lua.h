#pragma once

#include "framework/lua.h"

/** Creates a mesh handler and sets it as "current".
 *
 * \return Handle int Handle to the created mesh handler.
 * \ingroup LuaMeshHandler
 * \author Jan
 */
int MeshHandlerCreate(lua_State* L);

/**Exports the mesh to a wavefront.obj format.
 * \param FileName char Base name of the file to be used.
 * \param ExportByMaterial bool Default: False. Flag indicating whether to export
 *                      the extruder's surface mesh by material.
 * \ingroup LuaMeshHandler
 */
int MeshHandlerExportMeshToObj(lua_State* L);

/**Exports the mesh to vtu format.
 * \param FileName char Base name of the file to be used.
 * \ingroup LuaMeshHandler
 */
int MeshHandlerExportMeshToVTK(lua_State* L);

/**Exports the mesh to exodus format (.e extensions).
 * \param FileName char Base name of the file to be used.
 * \param suppress_nodesets bool Optional. Flag to suppress exporting nodesets.
 *                               Default = `false`.
 * \param suppress_sidesets bool Optional. Flag to suppress exporting sidesets.
 *                               Default = `false`.
 * \ingroup LuaMeshHandler
 */
int MeshHandlerExportMeshToExodus(lua_State* L);
