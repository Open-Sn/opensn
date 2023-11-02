#ifndef FV_DIFFUSION_LUA_UTILS_H
#define FV_DIFFUSION_LUA_UTILS_H

#include "chi_lua.h"
#include "../fv_diffusion_solver.h"

namespace fv_diffusion::fv_diffusion_lua_utils
{
int chiFVDiffusionSolverCreate(lua_State* L);
int chiFVDiffusionSetBCProperty(lua_State* L);

void RegisterLuaEntities(lua_State* L);
} // namespace fv_diffusion::fv_diffusion_lua_utils

#endif // FV_DIFFUSION_LUA_UTILS_H