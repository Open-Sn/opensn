#include "diffusion_lua.h"

#define LUA_FMACRO1(x) lua_register(L, #x, x)

namespace opensnlua::diffusion_solver
{

void
RegisterLuaEntities(lua_State* L)
{
  LUA_FMACRO1(chiDiffusionCreateSolver);
  LUA_FMACRO1(chiDiffusionInitialize);
  LUA_FMACRO1(chiDiffusionExecute);
  LUA_FMACRO1(chiDiffusionSetProperty);
}

} // namespace opensnlua::diffusion_solver
