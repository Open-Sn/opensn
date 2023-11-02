#ifndef CFEM_MG_DIFFUSION_LUA_UTILS_H
#define CFEM_MG_DIFFUSION_LUA_UTILS_H

#include "chi_lua.h"
#include "../mg_diffusion_solver.h"

namespace mg_diffusion::mgd_lua_utils
{
int chiCFEMMGDiffusionSolverCreate(lua_State* L);
int chiCFEMMGDiffusionSetBCProperty(lua_State* L);

void RegisterLuaEntities(lua_State* L);
} // namespace mg_diffusion::mgd_lua_utils

#endif // CFEM_MG_DIFFUSION_LUA_UTILS_H