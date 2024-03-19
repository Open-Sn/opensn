#include "framework/console/console.h"

#include "framework/lua.h"

#include "framework/utils/timer.h"

namespace opensnlua
{

/**Makes the program sleep for the specified time in milliseconds.
 * \param time int Time in milliseconds to sleep for.
 * */
int Sleep(lua_State* L);

RegisterLuaFunctionAsIs(Sleep);

int
Sleep(lua_State* L)
{
  const std::string fname = __FUNCTION__;
  const int num_args = lua_gettop(L);
  if (num_args != 1)
    LuaPostArgAmountError(fname, 1, num_args);

  const auto time_to_sleep = LuaArg<int64_t>(L, 1);

  opensn::Sleep(std::chrono::milliseconds(time_to_sleep));

  return 0;
}

} // namespace opensnlua
