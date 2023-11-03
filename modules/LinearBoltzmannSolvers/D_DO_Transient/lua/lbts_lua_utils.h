#pragma once

#include "chi_lua.h"

namespace lbs::lbts_lua_utils
{
int chiLBSCreateTransientSolver(lua_State* L);
int chiLBTSSetProperty(lua_State* L);
int chiLBTSGetProperty(lua_State* L);
int chiLBTSAdvanceTimeData(lua_State* L);

void RegisterLuaEntities(lua_State* L);
} // namespace lbs::lbts_lua_utils


