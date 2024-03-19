#include "framework/console/console.h"

#include "framework/lua.h"

#include "framework/runtime.h"

namespace opensnlua
{

/**
 * Gracefully exits OpenSn.
 * \param return_code int Return code, defaults to 0 (Success).
 */
int Exit(lua_State* L);

RegisterLuaFunctionAsIs(Exit);

int
Exit(lua_State* L)
{
  auto return_code = LuaArgOptional<int>(L, 1, EXIT_SUCCESS);
  opensn::Exit(return_code);
  return 0;
}

} // namespace opensnlua
