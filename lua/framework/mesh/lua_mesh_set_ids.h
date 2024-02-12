#pragma once

#include "framework/lua.h"

/**
 * Sets all cell-material id's to the supplied value.
 */
int MeshSetMatIDToAll(lua_State* L);

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
void SetMatIDFromLuaFunction(const std::string& lua_fname);

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
void SetBoundaryIDFromLuaFunction(const std::string& lua_fname);
