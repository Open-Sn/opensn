#include "framework/console/console.h"

#include "framework/lua.h"

#include "framework/runtime.h"

namespace opensnlua
{

/**Gracefully exits ChiTech.
 * \param return_code int Return code, defaults to 0 (Success).*/
int Exit(lua_State* L);

RegisterLuaFunctionAsIs(Exit);

int
Exit(lua_State* L)
{
  const int num_args = lua_gettop(L);

  int return_code = EXIT_SUCCESS;
  if (num_args >= 1)
  {
    LuaCheckIntegerValue(__FUNCTION__, L, 1);
    return_code = lua_tointeger(L, 1);
  }

  opensn::Exit(return_code);
  return 0;
}

} // namespace opensnlua
