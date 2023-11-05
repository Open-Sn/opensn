#pragma once

#include "opensn/framework/chi_lua.h"
#include "opensn/modules/MGDiffusion/mg_diffusion_solver.h"

namespace mg_diffusion::mgd_lua_utils
{
int chiCFEMMGDiffusionSolverCreate(lua_State* L);
int chiCFEMMGDiffusionSetBCProperty(lua_State* L);

void RegisterLuaEntities(lua_State* L);
} // namespace mg_diffusion::mgd_lua_utils
