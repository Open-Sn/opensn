#ifndef CFEM_DIFFUSION_LUA_UTILS_H
#define CFEM_DIFFUSION_LUA_UTILS_H

#include "chi_lua.h"
#include "../cfem_diffusion_solver.h"

namespace cfem_diffusion::cfem_diffusion_lua_utils
{
int chiCFEMDiffusionSolverCreate(lua_State* L);
int chiCFEMDiffusionSetBCProperty(lua_State* L);

void RegisterLuaEntities(lua_State* L);
} // namespace cfem_diffusion::cfem_diffusion_lua_utils

#endif // CFEM_DIFFUSION_LUA_UTILS_H