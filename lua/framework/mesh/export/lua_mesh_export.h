#pragma once

#include "framework/lua.h"

/**
 * Exports the mesh to a wavefront OBJ format.
 */
int MeshExportToObj(lua_State* L);

/**
 * Exports the mesh to vtu format.
 */
int MeshExportToVTK(lua_State* L);

/**
 * Exports the mesh to exodus format.
 */
int MeshExportToExodus(lua_State* L);
