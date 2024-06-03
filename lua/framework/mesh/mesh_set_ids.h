// SPDX-FileCopyrightText: 2024 The OpenSn Authors <https://open-sn.github.io/opensn/>
// SPDX-License-Identifier: MIT

#pragma once

#include "lua/framework/lua.h"

namespace opensnlua
{

/**
 * Sets all cell-material id's to the supplied value.
 */
int MeshSetUniformMaterialID(lua_State* L);

/**
 * Sets material id's using a lua function. The lua function is called with for each cell with 4
 * arguments, the cell's centroid x,y,z values and the cell's current material id.
 *
 * The lua function's prototype should be:
 * \code
 * function LuaFuncName(x,y,z,id)
 *   --stuff
 * end
 * \endcode
 */
int MeshSetMaterialIDFromLuaFunction(lua_State* L);

/**
 * Set specified material IDs using a LogicalVolume
 */
int MeshSetMaterialIDFromLogicalVolume(lua_State* L);

/**
 * Sets boundary id's using a lua function. The lua function is called for each boundary face
 * with 7 arguments, the face's centroid x,y,z values, the face's normal x,y,z values and the
 * face's current boundary id. The function must return a new_bndry_name (string).
 *
 * The lua function's prototype should be:
 * \code
 * function LuaFuncName(x,y,z,nx,ny,nz,id)
 * --stuff
 * end
 * \endcode
 */
int MeshSetBoundaryIDFromLuaFunction(lua_State* L);

/**
 * Set specified boundary IDs using a LogicalVolume
 */
int MeshSetBoundaryIDFromLogicalVolume(lua_State* L);

} // namespace opensnlua
