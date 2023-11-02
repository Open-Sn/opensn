//
// Created by Ragusa, Jean C on 12/18/22.
//
#include "ip_lua_utils.h"

#include "chi_runtime.h"

#define LUA_FMACRO1(x) lua_register(L, #x, x)
#define LUA_CMACRO1(x, y)                                                                          \
  lua_pushnumber(L, y);                                                                            \
  lua_setglobal(L, #x)

void
dfem_diffusion::dfem_diffusion_lua_utils::RegisterLuaEntities(lua_State* L)
{
  LUA_FMACRO1(chiDFEMDiffusionSolverCreate);
  LUA_FMACRO1(chiDFEMDiffusionSetBCProperty);

  LUA_CMACRO1(MAX_ITERATIONS, 1);
  LUA_CMACRO1(TOLERANCE, 2);
}