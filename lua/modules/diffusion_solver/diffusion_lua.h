#pragma once

#include "framework/chi_lua.h"

int chiDiffusionCreateSolver(lua_State* L);
int chiDiffusionInitialize(lua_State* L);
int chiDiffusionExecute(lua_State* L);
int chiDiffusionSetProperty(lua_State* L);

namespace diffusion_solver::lua_utils
{
void RegisterLuaEntities(lua_State* L);
} // namespace diffusion_solver::lua_utils
