#pragma once

#include "opensn/framework/chi_lua.h"
#include "opensn/modules/DFEMDiffusion/dfem_diffusion_solver.h"

int chiDFEMDiffusionSolverCreate(lua_State* L);
int chiDFEMDiffusionSetBCProperty(lua_State* L);

namespace dfem_diffusion
{
namespace dfem_diffusion_lua_utils
{
void RegisterLuaEntities(lua_State* L);
} // namespace dfem_diffusion_lua_utils
} // namespace dfem_diffusion
