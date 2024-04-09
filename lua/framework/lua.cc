#include "lua/framework/lua.h"

#include "framework/logging/log.h"

#include <string>
#include <sstream>
#include <map>

static int a = 15;

void
LuaPostArgAmountError(const std::string& func_name, int expected, int given)
{
  throw std::invalid_argument("Incorrect amount of arguments supplied in " + func_name +
                              " expected " + std::to_string(expected) + " arguments " + " but " +
                              std::to_string(given) + " provided");
}

void
LuaCheckNilValue(const std::string& func_name, lua_State* L, int arg)
{
  if (lua_isnil(L, arg))
  {
    throw std::invalid_argument("Nil -value supplied in " + func_name + " argument " +
                                std::to_string(arg));
  }
}

void
LuaCheckStringValue(const std::string& func_name, lua_State* L, int arg)
{
  if (not lua_isstring(L, arg))
  {
    throw std::invalid_argument("Non-string value supplied in " + func_name + " argument " +
                                std::to_string(arg));
  }
}

void
LuaCheckBoolValue(const std::string& func_name, lua_State* L, int arg)
{
  if (not lua_isboolean(L, arg))
  {
    throw std::invalid_argument("Non-boolean value supplied in " + func_name + " argument " +
                                std::to_string(arg));
  }
}

void
LuaCheckNumberValue(const std::string& func_name, lua_State* L, int arg)
{
  if (not lua_isnumber(L, arg))
  {
    throw std::invalid_argument("Non-number value supplied in " + func_name + " argument " +
                                std::to_string(arg));
  }
}

void
LuaCheckIntegerValue(const std::string& func_name, lua_State* L, int arg)
{
  if (not lua_isinteger(L, arg))
  {
    throw std::invalid_argument("Non-integer value supplied in " + func_name + " argument " +
                                std::to_string(arg));
  }
}

void
LuaCheckTableValue(const std::string& func_name, lua_State* L, int arg)
{
  if (not lua_istable(L, arg))
  {
    throw std::invalid_argument("Non-table value supplied in " + func_name + " argument " +
                                std::to_string(arg));
  }
}

std::string
LuaSourceInfo(lua_State* L, const char* func_name)
{
  lua_Debug err_info;
  lua_getstack(L, 1, &err_info);
  lua_getinfo(L, "nSl", &err_info);

  std::stringstream ret_str;
  ret_str << func_name << " " << err_info.source << " line " << err_info.currentline;

  return ret_str.str();
}
