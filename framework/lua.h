#pragma once

#ifdef OPENSN_WITH_LUA
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

namespace chi
{
class ParameterBlock;
}
namespace chi_lua
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
                                          chi::ParameterBlock& block,
                                          const std::string& key_str_name);

  /**This function operates on table keys recursively. It has a specific
   * behavior if it detects an array.*/
  static void RecursivelyParseTableKeys(lua_State* L, int t, chi::ParameterBlock& block);

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
  static chi::ParameterBlock ParseTable(lua_State* L, int table_stack_index);
};

/**\brief This recursive routine operates on a parameter block and passes
 * the parameters to the lua stack.
 *
 * If the `level` parameter is left as default then the zeroth level of
 * the parameter block will have its individual parameters exported as single
 * values, otherwise the block is exported as a table.
 */
void PushParameterBlock(lua_State* L, const chi::ParameterBlock& block, int level = 0);

chi::ParameterBlock StackItemToParameterBlock(lua_State* L, int index);
} // namespace chi_lua
#endif
