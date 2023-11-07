#ifdef OPENSN_WITH_LUA
#include "framework/chi_lua.h"

#include "framework/logging/chi_log.h"

#include <string>
#include <sstream>
#include <map>

static int a = 15;

/**Posts a generalized error message indicating that the
 * expected amount of arguments don't match the given amount.*/
void
LuaPostArgAmountError(const std::string& func_name, int expected, int given)
{
  throw std::invalid_argument("Incorrect amount of arguments supplied in " + func_name +
                              " expected " + std::to_string(expected) + " arguments " + " but " +
                              std::to_string(given) + " provided");
}

/**Checks if the lua variable at the stack location indicated by <arg>
 * is a nil value. Throws an error if it is.*/
void
LuaCheckNilValue(const std::string& func_name, lua_State* L, int arg)
{
  if (lua_isnil(L, arg))
  {
    throw std::invalid_argument("Nil -value supplied in " + func_name + " argument " +
                                std::to_string(arg));
  }
}

/**Checks if the lua variable at the stack location indicated by <arg>
 * is a string. Throws an error if it is not.*/
void
LuaCheckStringValue(const std::string& func_name, lua_State* L, int arg)
{
  if (not lua_isstring(L, arg))
  {
    throw std::invalid_argument("Non-string value supplied in " + func_name + " argument " +
                                std::to_string(arg));
  }
}

/**Checks if the lua variable at the stack location indicated by <arg>
 * is a boolean. Throws an error if it is not.*/
void
LuaCheckBoolValue(const std::string& func_name, lua_State* L, int arg)
{
  if (not lua_isboolean(L, arg))
  {
    throw std::invalid_argument("Non-boolean value supplied in " + func_name + " argument " +
                                std::to_string(arg));
  }
}

/**Checks if the lua variable at the stack location indicated by <arg>
 * is a number. Throws an error if it is not.*/
void
LuaCheckNumberValue(const std::string& func_name, lua_State* L, int arg)
{
  if (not lua_isnumber(L, arg))
  {
    throw std::invalid_argument("Non-number value supplied in " + func_name + " argument " +
                                std::to_string(arg));
  }
}

/**Checks if the lua variable at the stack location indicated by <arg>
 * is an integer. Throws an error if it is not.*/
void
LuaCheckIntegerValue(const std::string& func_name, lua_State* L, int arg)
{
  if (not lua_isinteger(L, arg))
  {
    throw std::invalid_argument("Non-integer value supplied in " + func_name + " argument " +
                                std::to_string(arg));
  }
}

/**Checks if the lua variable at the stack location indicated by <arg>
 * is an actual table. Throws an error if it is not.*/
void
LuaCheckTableValue(const std::string& func_name, lua_State* L, int arg)
{
  if (not lua_istable(L, arg))
  {
    throw std::invalid_argument("Non-table value supplied in " + func_name + " argument " +
                                std::to_string(arg));
  }
}

/**Gets information about an error state.*/
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

/**Performs the necessary checks and converts a lua-table, formatted
 * as a 1D array (i.e. numerical keys) to a std::vector. If the table
 * is not convertable, an error message is posted.*/
void
LuaPopulateVectorFrom1DArray(const std::string& func_name,
                             lua_State* L,
                             int table_arg_index,
                             std::vector<double>& vec)
{
  LuaCheckTableValue(func_name, L, table_arg_index);

  // Get the table as map
  std::map<int, double> vec_map;
  lua_pushnil(L);
  while (lua_next(L, table_arg_index) != 0)
  {
    if (not lua_isinteger(L, -2)) goto invalid_table;
    if (not lua_isnumber(L, -1)) goto invalid_table;
    int key = lua_tonumber(L, -2);
    double val = lua_tonumber(L, -1);

    vec_map[key] = val;

    lua_pop(L, 1);
  }

  // Populate vector
  {
    const size_t vec_size = vec_map.size();
    vec.clear();
    vec.resize(vec_size, 0.0);

    for (auto v_i : vec_map)
      vec.at(v_i.first - 1) = v_i.second;

    return;
  }

invalid_table:
  throw std::invalid_argument("Invalid table used in call to " + func_name);
}
#endif
