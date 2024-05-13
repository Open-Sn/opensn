// SPDX-FileCopyrightText: 2024 The OpenSn Authors <https://open-sn.github.io/opensn/>
// SPDX-License-Identifier: MIT

#pragma once

#include "framework/lua.h"

namespace opensnlua
{

/**
 * Exports the mesh to a wavefront OBJ format.
 */
int MeshExportToOBJ(lua_State* L);

/**
 * Exports the mesh to vtu format.
 */
int MeshExportToVTK(lua_State* L);

/**
 * Exports the mesh to ExodusII format.
 */
int MeshExportToExodusII(lua_State* L);

} // namespace opensnlua
