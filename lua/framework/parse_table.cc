// SPDX-FileCopyrightText: 2024 The OpenSn Authors <https://open-sn.github.io/opensn/>
// SPDX-License-Identifier: MIT

#include "lua/framework/lua.h"
#include "framework/logging/log_exceptions.h"
#include "framework/parameters/parameter_block.h"

namespace opensnlua
{

namespace
{

//  NOLINTBEGIN(misc-no-recursion)

/**
 * This function recursively processes table values. If the value is
 * a primitive type the recursion stops and the parameter block, which is
 * currently active, will be extended with a parameter of this primitive
 * type. If the value is another table, a new `Block`-type will be instantiated
 * and the table recursion will then operate on this new block.
 */
void
RecursivelyParseTableValues(lua_State* L,
                            opensn::ParameterBlock& block,
                            const std::string& key_str_name)
{
  switch (lua_type(L, -1))
  {
    case LUA_TNIL:
      throw std::logic_error(std::string(__PRETTY_FUNCTION__) +
                             ": Encountered nil value assigned to key " + key_str_name);
    case LUA_TBOOLEAN:
    {
      const bool bool_value = lua_toboolean(L, -1);
      block.AddParameter(key_str_name, bool_value);
      break;
    }
    case LUA_TNUMBER:
    {
      if (lua_isinteger(L, -1))
      {
        const int64_t number_value = lua_tointeger(L, -1);
        block.AddParameter(key_str_name, number_value);
      }
      else
      {
        const double number_value = lua_tonumber(L, -1);
        block.AddParameter(key_str_name, number_value);
      }

      break;
    }
    case LUA_TSTRING:
    {
      const std::string string_value = lua_tostring(L, -1);
      block.AddParameter(key_str_name, string_value);
      break;
    }
    case LUA_TTABLE:
    {
      opensn::ParameterBlock new_block(key_str_name);
      RecursivelyParseTableKeys(L, lua_gettop(L), new_block);
      block.AddParameter(new_block);
      break;
    }
    default:
      throw std::logic_error(std::string(__PRETTY_FUNCTION__) +
                             ": Encountered unsupported value type " +
                             lua_typename(L, lua_type(L, -2)) + " for key " + key_str_name);
  }
}
// NOLINTEND(misc-no-recursion)

} // namespace

// NOLINTBEGIN(misc-no-recursion)
void
RecursivelyParseTableKeys(lua_State* L, int t, opensn::ParameterBlock& block)
{
  bool number_key_encountered = false;
  bool string_key_encountered = false;

  int key_number_index = 0;

  lua_pushnil(L);             // first key
  while (lua_next(L, t) != 0) // pops the key, pushes next key and value
  {
    if (lua_type(L, -2) == LUA_TSTRING)
    {
      if (number_key_encountered)
        throw std::logic_error(std::string(__PRETTY_FUNCTION__) +
                               ": Encountered mixed key types (string and number)");

      string_key_encountered = true;
      const std::string key_str_name = lua_tostring(L, -2);
      RecursivelyParseTableValues(L, block, key_str_name);
    } // if key is string

    // If the key is a number then the following apply:
    // - This must be an array of items
    // - All the keys in the table must be numbers
    if (lua_type(L, -2) == LUA_TNUMBER)
    {
      if (string_key_encountered)
        throw std::logic_error(std::string(__PRETTY_FUNCTION__) +
                               ": Encountered mixed key types (string and number)");

      if (block.Type() != opensn::ParameterBlockType::ARRAY)
        block.ChangeToArray();

      number_key_encountered = true;
      const std::string key_str_name = std::to_string(key_number_index);
      RecursivelyParseTableValues(L, block, key_str_name);
      ++key_number_index;
    }

    lua_pop(L, 1);
  }
}
// NOLINTEND(misc-no-recursion)

//  NOLINTBEGIN(misc-no-recursion)
void
PushParameterBlock(lua_State* L, const opensn::ParameterBlock& block, int level)
{
  switch (block.Type())
  {
    case opensn::ParameterBlockType::BOOLEAN:
      LuaPush(L, block.GetValue<bool>());
      break;
    case opensn::ParameterBlockType::FLOAT:
      LuaPush(L, block.GetValue<double>());
      break;
    case opensn::ParameterBlockType::STRING:
      LuaPush(L, block.GetValue<std::string>());
      break;
    case opensn::ParameterBlockType::INTEGER:
      LuaPush(L, block.GetValue<lua_Integer>());
      break;
    case opensn::ParameterBlockType::ARRAY:
    {
      if (level > 0)
        lua_newtable(L);
      const size_t num_params = block.NumParameters();
      for (size_t k = 0; k < num_params; ++k)
      {
        if (level > 0)
          LuaPush(L, k + 1);
        PushParameterBlock(L, block.GetParam(k), level + 1);
        if (level > 0)
          lua_settable(L, -3);
      }
      break;
    }
    case opensn::ParameterBlockType::BLOCK:
    {
      if (level > 0)
        lua_newtable(L);
      const size_t num_params = block.NumParameters();
      for (size_t k = 0; k < num_params; ++k)
      {
        const auto& param = block.GetParam(k);
        if (level > 0)
          LuaPush(L, param.Name());
        PushParameterBlock(L, block.GetParam(k), level + 1);
        if (level > 0)
          lua_settable(L, -3);
      }
      break;
    }
    default:
      OpenSnLogicalError("Attempting to push unsupported ParameterBlockType to lua");
  }
}
//  NOLINTEND(misc-no-recursion)

} // namespace opensnlua
