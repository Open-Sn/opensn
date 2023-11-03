#pragma once

namespace chi
{
class InputParameters;
}

#include "chi_lua.h"
namespace chi_physics::lua_utils
{
int chiSolverCreate(lua_State* L);
int chiSolverInitialize(lua_State* L);
int chiSolverExecute(lua_State* L);
int chiSolverStep(lua_State* L);
int chiSolverAdvance(lua_State* L);
int chiSolverSetBasicOption(lua_State* L);
int chiSolverGetName(lua_State* L);
int chiSolverGetFieldFunctionList(lua_State* L);
int chiSolverGetInfo(lua_State* L);
int chiSolverSetProperties(lua_State* L);
} // namespace chi_physics::lua_utils
