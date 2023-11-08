#pragma once

#include "framework/chi_lua.h"
#include "modules/FVDiffusion/fv_diffusion_solver.h"

namespace fv_diffusion::fv_diffusion_lua_utils
{
int chiFVDiffusionSolverCreate(lua_State* L);
int chiFVDiffusionSetBCProperty(lua_State* L);

void RegisterLuaEntities(lua_State* L);
} // namespace fv_diffusion::fv_diffusion_lua_utils
