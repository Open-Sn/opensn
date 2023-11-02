#ifndef CHI_LUA_H
#define CHI_LUA_H

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

void LuaPostArgAmountError(const std::string& func_name, int expected, int given);
void LuaCheckNilValue(const std::string& func_name, lua_State* L, int arg);
void LuaCheckStringValue(const std::string& func_name, lua_State* L, int arg);
void LuaCheckBoolValue(const std::string& func_name, lua_State* L, int arg);
void LuaCheckNumberValue(const std::string& func_name, lua_State* L, int arg);
void LuaCheckIntegerValue(const std::string& func_name, lua_State* L, int arg);
void LuaCheckTableValue(const std::string& func_name, lua_State* L, int arg);
std::string LuaSourceInfo(lua_State* L, const char* func_name);
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
  static void RecursivelyParseTableValues(lua_State* L,
                                          chi::ParameterBlock& block,
                                          const std::string& key_str_name);

  static void RecursivelyParseTableKeys(lua_State* L, int t, chi::ParameterBlock& block);

public:
  /**\brief Parses a lua table into a hierarchical parameter block.*/
  static chi::ParameterBlock ParseTable(lua_State* L, int table_stack_index);
};

/**\brief This recursive routine operates on a parameter block and passes
 * the parameters to the lua stack.*/
void PushParameterBlock(lua_State* L, const chi::ParameterBlock& block, int level = 0);

chi::ParameterBlock StackItemToParameterBlock(lua_State* L, int index);
} // namespace chi_lua

#endif