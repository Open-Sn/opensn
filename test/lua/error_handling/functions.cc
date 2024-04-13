#include "framework/runtime.h"
#include "lua/framework/lua.h"
#include "lua/framework/console/console.h"
#include "framework/logging/log.h"

using namespace opensnlua;

namespace unit_tests
{

int
FunctionWithCheck(lua_State* L)
{
  const std::string fname = "FunctionWithCheck";
  LuaCheckArgs<std::string, int>(L, fname);
  return LuaReturn(L);
}

int
FunctionWithoutCheck(lua_State* L)
{
  auto str = LuaArg<std::string>(L, 1);
  auto i = LuaArg<int>(L, 2);
  opensn::log.Log() << "str = " << str << ", i = " << i;
  return LuaReturn(L);
}

RegisterLuaFunctionNamespace(FunctionWithCheck, unit_tests, FunctionWithCheck);
RegisterLuaFunctionNamespace(FunctionWithoutCheck, unit_tests, FunctionWithoutCheck);

} //  namespace unit_tests
