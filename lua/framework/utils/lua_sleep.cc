#include "framework/console/console.h"

#include "framework/lua.h"

#include "framework/utils/timer.h"

namespace chi::lua_utils
{

/**Makes the program sleep for the specified time in milliseconds.
 * \param time int Time in milliseconds to sleep for.
 * */
int chiSleep(lua_State* L);

RegisterLuaFunctionAsIs(chiSleep);

int
chiSleep(lua_State* L)
{
  const std::string fname = __FUNCTION__;
  const int num_args = lua_gettop(L);
  if (num_args != 1) LuaPostArgAmountError(fname, 1, num_args);

  LuaCheckIntegerValue(fname, L, 1);
  const int64_t time_to_sleep = lua_tointeger(L, 1);

  chi::Sleep(std::chrono::milliseconds(time_to_sleep));

  return 0;
}

} // namespace chi::lua_utils
