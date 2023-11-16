#ifdef OPENSN_WITH_LUA
#include "framework/lua.h"

#include "framework/runtime.h"
#include "framework/logging/log_exceptions.h"

#include "framework/data_types/data_types.h"
#include "framework/parameters/parameter_block.h"

#define ExceptionLuaNilValue                                                                       \
  throw std::logic_error(std::string(__PRETTY_FUNCTION__) +                                        \
                         ": Encountered nil value assigned to key " + key_str_name)

#define ExceptionLuaUnsupportedValue                                                               \
  throw std::logic_error(std::string(__PRETTY_FUNCTION__) +                                        \
                         ": Encountered unsupported value type " +                                 \
                         lua_typename(L, lua_type(L, -2)) + " for key " + key_str_name)

#define ExceptionMixStringNumberKeys                                                               \
  throw std::logic_error(std::string(__PRETTY_FUNCTION__) +                                        \
                         ": Encountered mixed key types (string and number)")

namespace opensn
{

typedef ParameterBlock ParamBlock;

//  NOLINTBEGIN(misc-no-recursion)
void
TableParserAsParameterBlock::RecursivelyParseTableValues(lua_State* L,
                                                         ParamBlock& block,
                                                         const std::string& key_str_name)
{
  switch (lua_type(L, -1))
  {
    case LUA_TNIL:
      ExceptionLuaNilValue;
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
      ParameterBlock new_block(key_str_name);
      RecursivelyParseTableKeys(L, lua_gettop(L), new_block);
      block.AddParameter(new_block);
      break;
    }
    default:
      ExceptionLuaUnsupportedValue;
  } // switch on value types
}
// NOLINTEND(misc-no-recursion)

//  NOLINTBEGIN(misc-no-recursion)
void
TableParserAsParameterBlock::RecursivelyParseTableKeys(lua_State* L, int t, ParameterBlock& block)
{
  bool number_key_encountered = false;
  bool string_key_encountered = false;

  int key_number_index = 0;

  lua_pushnil(L);             // first key
  while (lua_next(L, t) != 0) // pops the key, pushes next key and value
  {
    if (lua_type(L, -2) == LUA_TSTRING)
    {
      if (number_key_encountered) ExceptionMixStringNumberKeys;

      string_key_encountered = true;
      const std::string key_str_name = lua_tostring(L, -2);
      RecursivelyParseTableValues(L, block, key_str_name);
    } // if key is string

    // If the key is a number then the following apply:
    // - This must be an array of items
    // - All the keys in the table must be numbers
    if (lua_type(L, -2) == LUA_TNUMBER)
    {
      if (string_key_encountered) ExceptionMixStringNumberKeys;

      if (block.Type() != ParameterBlockType::ARRAY) block.ChangeToArray();

      number_key_encountered = true;
      const std::string key_str_name = std::to_string(key_number_index);
      RecursivelyParseTableValues(L, block, key_str_name);
      ++key_number_index;
    }

    lua_pop(L, 1);
  }
}
// NOLINTEND(misc-no-recursion)

ParameterBlock
TableParserAsParameterBlock::ParseTable(lua_State* L, int table_stack_index)
{
  ParamBlock param_block;

  RecursivelyParseTableKeys(L, table_stack_index, param_block);

  return param_block;
}

//  NOLINTBEGIN(misc-no-recursion)
void
PushParameterBlock(lua_State* L, const ParameterBlock& block, int level)
{
  switch (block.Type())
  {
    case ParameterBlockType::BOOLEAN:
      lua_pushboolean(L, block.GetValue<bool>());
      break;
    case ParameterBlockType::FLOAT:
      lua_pushnumber(L, block.GetValue<double>());
      break;
    case ParameterBlockType::STRING:
      lua_pushstring(L, block.GetValue<std::string>().c_str());
      break;
    case ParameterBlockType::INTEGER:
      lua_pushinteger(L, block.GetValue<lua_Integer>());
      break;
    case ParameterBlockType::ARRAY:
    {
      if (level > 0) lua_newtable(L);
      const size_t num_params = block.NumParameters();
      for (size_t k = 0; k < num_params; ++k)
      {
        if (level > 0) lua_pushinteger(L, static_cast<lua_Integer>(k) + 1);
        PushParameterBlock(L, block.GetParam(k), level + 1);
        if (level > 0) lua_settable(L, -3);
      }
      break;
    }
    case ParameterBlockType::BLOCK:
    {
      if (level > 0) lua_newtable(L);
      const size_t num_params = block.NumParameters();
      for (size_t k = 0; k < num_params; ++k)
      {
        const auto& param = block.GetParam(k);
        if (level > 0) lua_pushstring(L, param.Name().c_str());
        PushParameterBlock(L, block.GetParam(k), level + 1);
        if (level > 0) lua_settable(L, -3);
      }
      break;
    }
    default:
      ChiLogicalError("Attempting to push unsupport ParameterBlockType to lua");
  }
}
//  NOLINTEND(misc-no-recursion)

ParameterBlock
StackItemToParameterBlock(lua_State* L, int index)
{
  switch (lua_type(L, index))
  {
    case LUA_TNIL:
      return ParamBlock{};
    case LUA_TBOOLEAN:
      return ParamBlock("", lua_toboolean(L, index));
    case LUA_TNUMBER:
    {
      if (lua_isinteger(L, index)) return ParamBlock("", lua_tointeger(L, index));
      else
        return ParamBlock("", lua_tonumber(L, index));
    }
    case LUA_TSTRING:
    {
      const std::string value = lua_tostring(L, index);
      return ParamBlock("", value);
    }
    case LUA_TTABLE:
    {
      auto paramblock = TableParserAsParameterBlock::ParseTable(L, index);
      return paramblock;
    }
    default:
      ChiLogicalError("Unhandled Lua type.");
  }
}

} // namespace opensn
#endif
