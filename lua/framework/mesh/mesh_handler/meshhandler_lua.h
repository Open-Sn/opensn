#pragma once

#include "framework/chi_lua.h"

/** Creates a mesh handler and sets it as "current".
 *
 * \return Handle int Handle to the created mesh handler.
 * \ingroup LuaMeshHandler
 * \author Jan
 */
int chiMeshHandlerCreate(lua_State* L);

/** Sets the given mesh handler as "current".
 *
 * \param HandlerHandler int Handle to the mesh handler previously created
 *        with a call to chiMeshHandlerCreate.
 *
 * \ingroup LuaMeshHandler
 * \author Jan
 */
int chiMeshHandlerSetCurrent(lua_State* L);

/**Exports the mesh to a wavefront.obj format.
 * \param FileName char Base name of the file to be used.
 * \param ExportByMaterial bool Default: False. Flag indicating whether to export
 *                      the extruder's surface mesh by material.
 * \ingroup LuaMeshHandler
 */
int chiMeshHandlerExportMeshToObj(lua_State* L);

/**Exports the mesh to vtu format.
 * \param FileName char Base name of the file to be used.
 * \ingroup LuaMeshHandler
 */
int chiMeshHandlerExportMeshToVTK(lua_State* L);

/**Exports the mesh to exodus format (.e extensions).
 * \param FileName char Base name of the file to be used.
 * \param suppress_nodesets bool Optional. Flag to suppress exporting nodesets.
 *                               Default = `false`.
 * \param suppress_sidesets bool Optional. Flag to suppress exporting sidesets.
 *                               Default = `false`.
 * \ingroup LuaMeshHandler
 */
int chiMeshHandlerExportMeshToExodus(lua_State* L);
