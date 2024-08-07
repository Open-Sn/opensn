// SPDX-FileCopyrightText: 2024 The OpenSn Authors <https://open-sn.github.io/opensn/>
// SPDX-License-Identifier: MIT

#pragma once

#include "lua/framework/lua.h"

namespace opensnlua
{
int AdjointSolverCreate(lua_State* L);
int AdjointSolverAddResponseFunction(lua_State* L);
int AdjointSolverMakeExpRepFromP1Moments(lua_State* L);
int AdjointSolverExportImportanceMapBinary(lua_State* L);
} // namespace opensnlua
