#include "diffusion_lua.h"

#define LUA_FMACRO1(x) lua_register(L, #x, x)

namespace opensnlua::diffusion_solver
{

void
RegisterLuaEntities(lua_State* L)
{
  LUA_FMACRO1(DiffusionCreateSolver);
  LUA_FMACRO1(DiffusionInitialize);
  LUA_FMACRO1(DiffusionExecute);
  LUA_FMACRO1(DiffusionSetProperty);
}

} // namespace opensnlua::diffusion_solver
