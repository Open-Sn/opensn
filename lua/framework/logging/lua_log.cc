#include "lua_log.h"
#include "framework/logging/log.h"
#include "framework/console/console.h"
#include "framework/lua.h"
#include "framework/runtime.h"

using namespace opensn;

namespace opensnlua
{
RegisterLuaFunctionNamespace(LogSetVerbosity, log, SetVerbosity);
RegisterLuaFunctionNamespace(LogLog, log, Log);

RegisterLuaConstantAsIs(LOG_0, Varying(1));
RegisterLuaConstantAsIs(LOG_0WARNING, Varying(2));
RegisterLuaConstantAsIs(LOG_0ERROR, Varying(3));
RegisterLuaConstantAsIs(LOG_0VERBOSE_0, Varying(4));
RegisterLuaConstantAsIs(LOG_0VERBOSE_1, Varying(5));
RegisterLuaConstantAsIs(LOG_0VERBOSE_2, Varying(6));
RegisterLuaConstantAsIs(LOG_ALL, Varying(7));
RegisterLuaConstantAsIs(LOG_ALLWARNING, Varying(8));
RegisterLuaConstantAsIs(LOG_ALLERROR, Varying(9));
RegisterLuaConstantAsIs(LOG_ALLVERBOSE_0, Varying(10));
RegisterLuaConstantAsIs(LOG_ALLVERBOSE_1, Varying(11));
RegisterLuaConstantAsIs(LOG_ALLVERBOSE_2, Varying(12));

int
LogSetVerbosity(lua_State* L)
{
  int num_args = lua_gettop(L);

  if (num_args == 0)
  {
    return 0;
  }
  else
  {
    int level = lua_tonumber(L, 1);
    if (level <= 2)
    {
      opensn::log.SetVerbosity(level);
    }
  }
  return 0;
}

int
LogLog(lua_State* L)
{
  int num_args = lua_gettop(L);

  if (num_args != 2)
    LuaPostArgAmountError("Log", 2, num_args);

  int mode = lua_tonumber(L, 1);
  const char* message = lua_tostring(L, 2);

  opensn::log.Log(static_cast<opensn::Logger::LOG_LVL>(mode)) << message << std::endl;

  return 0;
}

} // namespace opensnlua
