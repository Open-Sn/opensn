#include "framework/console/chi_console.h"

#include "framework/chi_lua.h"

#include "framework/chi_runtime.h"

namespace chi::lua_utils
{

/**Gracefully exits ChiTech.
 * \param return_code int Return code, defaults to 0 (Success).*/
int chiExit(lua_State* L);

RegisterLuaFunctionAsIs(chiExit);

int
chiExit(lua_State* L)
{
  const int num_args = lua_gettop(L);

  int return_code = EXIT_SUCCESS;
  if (num_args >= 1)
  {
    LuaCheckIntegerValue(__FUNCTION__, L, 1);
    return_code = lua_tointeger(L, 1);
  }

  Chi::Exit(return_code);
  return 0;
}

} // namespace chi::lua_utils
