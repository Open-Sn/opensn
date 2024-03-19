#pragma once

extern "C"
{
#include <lua.h>
#include <lualib.h>
#include <lauxlib.h>
}

#include <typeinfo>
#include <string>
#include <vector>
#include <memory>
#include "framework/parameters/parameter_block.h"

/**Posts a generalized error message indicating that the
 * expected amount of arguments don't match the given amount.*/
void LuaPostArgAmountError(const std::string& func_name, int expected, int given);
/**Checks if the lua variable at the stack location indicated by <arg>
 * is a nil value. Throws an error if it is.*/
void LuaCheckNilValue(const std::string& func_name, lua_State* L, int arg);
/**Checks if the lua variable at the stack location indicated by <arg>
 * is a string. Throws an error if it is not.*/
void LuaCheckStringValue(const std::string& func_name, lua_State* L, int arg);
/**Checks if the lua variable at the stack location indicated by <arg>
 * is a boolean. Throws an error if it is not.*/
void LuaCheckBoolValue(const std::string& func_name, lua_State* L, int arg);
/**Checks if the lua variable at the stack location indicated by <arg>
 * is a number. Throws an error if it is not.*/
void LuaCheckNumberValue(const std::string& func_name, lua_State* L, int arg);
/**Checks if the lua variable at the stack location indicated by <arg>
 * is an integer. Throws an error if it is not.*/
void LuaCheckIntegerValue(const std::string& func_name, lua_State* L, int arg);
/**Checks if the lua variable at the stack location indicated by <arg>
 * is an actual table. Throws an error if it is not.*/
void LuaCheckTableValue(const std::string& func_name, lua_State* L, int arg);
/**Gets information about an error state.*/
std::string LuaSourceInfo(lua_State* L, const char* func_name);
/**Performs the necessary checks and converts a lua-table, formatted
 * as a 1D array (i.e. numerical keys) to a std::vector. If the table
 * is not convertable, an error message is posted.*/
void LuaPopulateVectorFrom1DArray(const std::string& func_name,
                                  lua_State* L,
                                  int table_arg_index,
                                  std::vector<double>& vec);

namespace opensnlua
{

/**This static object is used to parse lua tables into parameter blocks.*/
class TableParserAsParameterBlock
{
private:
  /**This function recursively processes table values. If the value is
   * a primitive type the recursion stops and the parameter block, which is
   * currently active, will be extended with a parameter of this primitive
   * type. If the value is another table, a new `Block`-type will be instantiated
   * and the table recursion will then operate on this new block.*/
  static void RecursivelyParseTableValues(lua_State* L,
                                          opensn::ParameterBlock& block,
                                          const std::string& key_str_name);

  /**This function operates on table keys recursively. It has a specific
   * behavior if it detects an array.*/
  static void RecursivelyParseTableKeys(lua_State* L, int t, opensn::ParameterBlock& block);

public:
  /**\brief Parses a lua table into a hierarchical parameter block.
   *
   * This is the root command for parsing a table as a parameter block.
   * * Example table:
   * \code
   * block =
   * {
   *   enabled = true,
   *   it_method = "gmres",
   *   nl_abs_tol = 1.0e-12,
   *   nl_max_its = 33,
   *   sub1 =
   *   {
   *     ax_method = 2,
   *     l_abs_tol = 1.0e-2
   *   },
   *   sub2 =
   *   {
   *     ax_method = 3,
   *     l_abs_tol = 1.0e-3,
   *     blocks = {99, 98, 97},
   *     cblocks = {{1,2,3},{4,5,6},{7,8,9}}
   *   }
   * }
   *
   * chiUnitTests_Test_paramblock(--[[verbose=]]true, block)
   * \endcode
   */
  static opensn::ParameterBlock ParseTable(lua_State* L, int table_stack_index);
};

/**\brief This recursive routine operates on a parameter block and passes
 * the parameters to the lua stack.
 *
 * If the `level` parameter is left as default then the zeroth level of
 * the parameter block will have its individual parameters exported as single
 * values, otherwise the block is exported as a table.
 */
void PushParameterBlock(lua_State* L, const opensn::ParameterBlock& block, int level = 0);

opensn::ParameterBlock StackItemToParameterBlock(lua_State* L, int index);

template <typename T>
inline T
LuaArgAsType(lua_State* L, int index)
{
  throw std::invalid_argument("Unsupported type: " + std::string(typeid(T).name()));
}

template <>
inline bool
LuaArgAsType(lua_State* L, int index)
{
  if (lua_isboolean(L, index))
    return lua_toboolean(L, index) == 1;
  else
    throw std::invalid_argument("Expected boolean value as " + std::to_string(index) +
                                ". argument");
}

template <>
inline int
LuaArgAsType(lua_State* L, int index)
{
  if (lua_isinteger(L, index))
    return lua_tointeger(L, index);
  else
    throw std::invalid_argument("Expected int value as " + std::to_string(index) + ". argument");
}

template <>
inline int64_t
LuaArgAsType(lua_State* L, int index)
{
  if (lua_isinteger(L, index))
    return lua_tointeger(L, index);
  else
    throw std::invalid_argument("Expected int value as " + std::to_string(index) + ". argument");
}

template <>
inline double
LuaArgAsType(lua_State* L, int index)
{
  if (lua_isnumber(L, index))
    return lua_tonumber(L, index);
  else
    throw std::invalid_argument("Expected number value as " + std::to_string(index) + ". argument");
}

template <>
inline std::size_t
LuaArgAsType(lua_State* L, int index)
{
  if (lua_isinteger(L, index))
    return lua_tointeger(L, index);
  else
    throw std::invalid_argument("Expected std::size_t value as " + std::to_string(index) +
                                ". argument");
}

template <>
inline const char*
LuaArgAsType(lua_State* L, int index)
{
  if (lua_isstring(L, index))
    return lua_tostring(L, index);
  else
    throw std::invalid_argument("Expected string value as " + std::to_string(index) + ". argument");
}

template <>
inline std::string
LuaArgAsType(lua_State* L, int index)
{
  if (lua_isstring(L, index))
    return std::string(lua_tostring(L, index));
  else
    throw std::invalid_argument("Expected string value as " + std::to_string(index) + ". argument");
}

/**
 * Get argument with specified index from Lua stack
 *
 * \param L Lua stack
 * \param index Argument index from the lua stack (1-based)
 */
template <typename T>
inline T
LuaArg(lua_State* L, int index)
{
  auto num_args = lua_gettop(L);
  if (index > num_args)
    throw std::invalid_argument("Invalid argument index " + std::to_string(index) +
                                ". Supplied only " + std::to_string(num_args) + " arguments.");
  if (lua_isnil(L, index))
    throw std::invalid_argument("Unexpected value supplied as " + std::to_string(index) +
                                ". argument.");

  return LuaArgAsType<T>(L, index);
}

/**
 * Get optional argument with specified index from Lua stack
 *
 * \param L Lua stack
 * \param index Argument index from the lua stack (1-based)
 * \param default_value Default value to use if argument was not specified
 */
template <typename T>
inline T
LuaArgOptional(lua_State* L, int index, T default_value)
{
  const int num_args = lua_gettop(L);
  if (index <= num_args)
    return LuaArg<T>(L, index);
  else
    return default_value;
}

} // namespace opensnlua
