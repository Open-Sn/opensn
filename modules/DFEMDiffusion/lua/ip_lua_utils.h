#ifndef DFEM_DIFFUSION_LUA_UTILS_H
#define DFEM_DIFFUSION_LUA_UTILS_H

#include "chi_lua.h"
#include "../dfem_diffusion_solver.h"

int chiDFEMDiffusionSolverCreate(lua_State* L);
int chiDFEMDiffusionSetBCProperty(lua_State* L);

namespace dfem_diffusion
{
namespace dfem_diffusion_lua_utils
{
void RegisterLuaEntities(lua_State* L);
} // namespace dfem_diffusion_lua_utils
} // namespace dfem_diffusion

#endif // DFEM_DIFFUSION_LUA_UTILS_H