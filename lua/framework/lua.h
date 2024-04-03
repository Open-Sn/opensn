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
#include <type_traits>
#include "framework/parameters/parameter_block.h"
#include "framework/math/math.h"

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

namespace opensnlua
{

/**
 * This function operates on table keys recursively. It has a specific
 * behavior if it detects an array.
 */

void RecursivelyParseTableKeys(lua_State* L, int t, opensn::ParameterBlock& block);

/**
 * This recursive routine operates on a parameter block and passes the parameters to the lua stack.
 *
 * If the `level` parameter is left as default then the zeroth level of
 * the parameter block will have its individual parameters exported as single
 * values, otherwise the block is exported as a table.
 */
void PushParameterBlock(lua_State* L, const opensn::ParameterBlock& block, int level = 0);

//

template <typename>
struct is_std_vector : std::false_type
{
};

template <typename T, typename A>
struct is_std_vector<std::vector<T, A>> : std::true_type
{
};

template <typename>
struct is_std_map : std::false_type
{
};

template <class KEY, class T, class C, class A>
struct is_std_map<std::map<KEY, T, C, A>> : std::true_type
{
};

// Push values on stack

template <typename T>
inline void
LuaPush(lua_State* L, T value)
{
  throw std::invalid_argument("Unsupported type: " + std::string(typeid(T).name()));
}

template <>
inline void
LuaPush(lua_State* L, char* value)
{
  lua_pushstring(L, value);
}

template <>
inline void
LuaPush(lua_State* L, const char* value)
{
  lua_pushstring(L, value);
}

template <>
inline void
LuaPush(lua_State* L, char value)
{
  lua_pushinteger(L, static_cast<lua_Integer>(value));
}

template <>
inline void
LuaPush(lua_State* L, int value)
{
  lua_pushinteger(L, static_cast<lua_Integer>(value));
}

template <>
inline void
LuaPush(lua_State* L, long value)
{
  lua_pushinteger(L, static_cast<lua_Integer>(value));
}

template <>
inline void
LuaPush(lua_State* L, long long value)
{
  lua_pushinteger(L, static_cast<lua_Integer>(value));
}

template <>
inline void
LuaPush(lua_State* L, unsigned char value)
{
  lua_pushinteger(L, static_cast<lua_Integer>(value));
}

template <>
inline void
LuaPush(lua_State* L, unsigned int value)
{
  lua_pushinteger(L, static_cast<lua_Integer>(value));
}

template <>
inline void
LuaPush(lua_State* L, unsigned long value)
{
  lua_pushinteger(L, static_cast<lua_Integer>(value));
}

template <>
inline void
LuaPush(lua_State* L, unsigned long long value)
{
  lua_pushinteger(L, static_cast<lua_Integer>(value));
}

template <>
inline void
LuaPush(lua_State* L, float value)
{
  lua_pushnumber(L, static_cast<lua_Number>(value));
}

template <>
inline void
LuaPush(lua_State* L, double value)
{
  lua_pushnumber(L, static_cast<lua_Number>(value));
}

template <>
inline void
LuaPush(lua_State* L, bool value)
{
  lua_pushboolean(L, static_cast<lua_Integer>(value));
}

inline void
LuaPush(lua_State* L, const std::string& value)
{
  lua_pushstring(L, value.c_str());
}

inline void
LuaPush(lua_State* L, const opensn::Vector3& value)
{
  lua_newtable(L);
  LuaPush(L, "x");
  LuaPush(L, value.x);
  lua_settable(L, -3);
  LuaPush(L, "y");
  LuaPush(L, value.y);
  lua_settable(L, -3);
  LuaPush(L, "z");
  LuaPush(L, value.z);
  lua_settable(L, -3);
}

inline void
LuaPush(lua_State* L, const opensn::ParameterBlock& value)
{
  PushParameterBlock(L, value);
}

template <typename T1, typename T2>
inline void
LuaPush(lua_State* L, const std::pair<T1, T2>& value)
{
  lua_newtable(L);

  LuaPush(L, "first");
  LuaPush(L, value.first);
  lua_settable(L, -3);

  LuaPush(L, "second");
  LuaPush(L, value.second);
  lua_settable(L, -3);
}

template <typename T, typename A>
inline void
LuaPush(lua_State* L, const std::vector<T, A>& value)
{
  lua_newtable(L);
  for (size_t i = 0; i < value.size(); ++i)
  {
    LuaPush(L, i + 1);
    LuaPush(L, value[i]);
    lua_settable(L, -3);
  }
}

template <class KEY, class T, class C, class A>
inline void
LuaPush(lua_State* L, const std::map<KEY, T, C, A>& value)
{
  lua_newtable(L);
  for (auto& [k, v] : value)
  {
    LuaPush(L, k);
    LuaPush(L, v);
    lua_settable(L, -3);
  }
}

template <typename KEY, typename VAL>
inline void
LuaPushTableKey(lua_State* L, const KEY& key, const VAL& value)
{
  LuaPush(L, key);
  LuaPush(L, value);
  lua_settable(L, -3);
}

namespace
{

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

template <>
inline opensn::Vector3
LuaArgAsType(lua_State* L, int index)
{
  if (lua_istable(L, index))
  {
    opensn::Vector3 vec;

    LuaPush(L, "x");
    lua_gettable(L, index);
    vec.x = LuaArgAsType<double>(L, -1);
    lua_pop(L, 1);

    LuaPush(L, "y");
    lua_gettable(L, index);
    if (not lua_isnil(L, -1))
      vec.y = LuaArgAsType<double>(L, -1);
    lua_pop(L, 1);

    LuaPush(L, "z");
    lua_gettable(L, index);
    if (not lua_isnil(L, -1))
      vec.z = LuaArgAsType<double>(L, -1);
    lua_pop(L, 1);

    return vec;
  }
  else
    throw std::invalid_argument("Expected table value as " + std::to_string(index) + ". argument");
}

template <>
inline opensn::ParameterBlock
LuaArgAsType(lua_State* L, int index)
{
  if (lua_istable(L, index))
  {
    opensn::ParameterBlock param_block;
    RecursivelyParseTableKeys(L, index, param_block);
    return param_block;
  }
  else
    throw std::invalid_argument("Expected table value as " + std::to_string(index) + ". argument");
}

template <typename T, typename A>
inline std::vector<T, A>
LuaArgAsType(lua_State* L, int index)
{
  if (lua_istable(L, index))
  {
    std::vector<T, A> values;
    const size_t sz = lua_rawlen(L, index);
    values.resize(sz);
    for (size_t i = 0; i < sz; i++)
    {
      LuaPush(L, i + 1);
      lua_gettable(L, index);
      values[i] = LuaArgAsType<T>(L, -1);
      lua_pop(L, 1);
    }
    return values;
  }
  else
    throw std::invalid_argument("Expected table value as " + std::to_string(index) + ". argument");
}

template <class KEY, class T, class C, class A>
inline std::map<KEY, T, C, A>
LuaArgAsType(lua_State* L, int index)
{
  if (lua_istable(L, index))
  {
    std::map<KEY, T, C, A> values;
    lua_pushnil(L);
    while (lua_next(L, index) != 0)
    {
      auto k = LuaArgAsType<KEY>(L, -2);
      auto val = LuaArgAsType<T>(L, -1);
      values.insert(std::pair<KEY, T>(k, val));
      lua_pop(L, 1);
    }
    return values;
  }
  else
    throw std::invalid_argument("Expected table value as " + std::to_string(index) + ". argument");
}

template <typename T>
inline T
LuaArgImpl(lua_State* L, int index)
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

template <typename T, typename A>
inline std::vector<T, A>
LuaArgStdVectorImpl(lua_State* L, int index)
{
  auto num_args = lua_gettop(L);
  if (index > num_args)
    throw std::invalid_argument("Invalid argument index " + std::to_string(index) +
                                ". Supplied only " + std::to_string(num_args) + " arguments.");
  if (lua_isnil(L, index))
    throw std::invalid_argument("Unexpected value supplied as " + std::to_string(index) +
                                ". argument.");

  return LuaArgAsType<T, A>(L, index);
}

template <class KEY, class T, class C, class A>
inline std::map<KEY, T, C, A>
LuaArgStdMapImpl(lua_State* L, int index)
{
  auto num_args = lua_gettop(L);
  if (index > num_args)
    throw std::invalid_argument("Invalid argument index " + std::to_string(index) +
                                ". Supplied only " + std::to_string(num_args) + " arguments.");
  if (lua_isnil(L, index))
    throw std::invalid_argument("Unexpected value supplied as " + std::to_string(index) +
                                ". argument.");

  return LuaArgAsType<KEY, T, C, A>(L, index);
}

} // namespace

/**
 * Get argument with specified index from Lua stack
 *
 * \param L Lua stack
 * \param index Argument index from the lua stack (1-based)
 */
template <typename TRET>
inline TRET
LuaArg(lua_State* L, int index)
{
  if constexpr (is_std_vector<TRET>::value)
  {
    using T = typename TRET::value_type;
    return LuaArgStdVectorImpl<T, std::allocator<T>>(L, index);
  }
  else if constexpr (is_std_map<TRET>::value)
  {
    using KEY = typename TRET::key_type;
    using T = typename TRET::mapped_type;
    return LuaArgStdMapImpl<KEY, T, std::less<KEY>, std::allocator<std::pair<const KEY, T>>>(L,
                                                                                             index);
  }
  else
    return LuaArgImpl<TRET>(L, index);
}

namespace
{

template <typename T>
inline T
LuaArgOptionalImpl(lua_State* L, int index, T default_value)
{
  const int num_args = lua_gettop(L);
  if (index <= num_args)
    return LuaArgImpl<T>(L, index);
  else
    return default_value;
}

template <typename T, typename A>
inline std::vector<T, A>
LuaArgOptionalStdVectorImpl(lua_State* L, int index, const std::vector<T>& default_value)
{
  const int num_args = lua_gettop(L);
  if (index <= num_args)
    return LuaArgStdVectorImpl<T, A>(L, index);
  else
    return default_value;
}

} // namespace

/**
 * Get optional argument with specified index from Lua stack
 *
 * \param L Lua stack
 * \param index Argument index from the lua stack (1-based)
 * \param default_value Default value to use if argument was not specified
 */
template <typename TRET>
inline TRET
LuaArgOptional(lua_State* L, int index, TRET default_value)
{
  if constexpr (is_std_vector<TRET>::value)
  {
    using T = typename TRET::value_type;
    return LuaArgOptionalStdVectorImpl<T, std::allocator<T>>(L, index, default_value);
  }
  else
    return LuaArgOptionalImpl<TRET>(L, index, default_value);
}

// Set global values

inline void
LuaSetGlobal(lua_State* L, const std::string& name, const std::string& value)
{
  lua_pushstring(L, value.c_str());
  lua_setglobal(L, name.c_str());
}

inline void
LuaSetGlobal(lua_State* L, const std::string& name, const char* value)
{
  lua_pushstring(L, value);
  lua_setglobal(L, name.c_str());
}

inline void
LuaSetGlobal(lua_State* L, const std::string& name, int value)
{
  lua_pushinteger(L, value);
  lua_setglobal(L, name.c_str());
}

inline void
LuaSetGlobal(lua_State* L, const std::string& name, double value)
{
  lua_pushnumber(L, value);
  lua_setglobal(L, name.c_str());
}

// Support methods for pushing arguments on a Lua stack

namespace
{

inline void
LuaPushArgs(lua_State* L)
{
}

template <typename T, typename... ARGS>
inline void
LuaPushArgs(lua_State* L, T&& val, ARGS&&... args)
{
  LuaPush(L, val);
  LuaPushArgs(L, std::forward<ARGS>(args)...);
}

template <typename TRET, typename... ARGS>
inline TRET
LuaCallImpl(lua_State* L, const std::string& fn_name, ARGS... args)
{
  TRET ret_val;
  lua_getglobal(L, fn_name.c_str());
  if (not lua_isfunction(L, -1))
    throw std::logic_error("Lua function '" + fn_name + "' does not exist.");

  LuaPushArgs(L, std::forward<ARGS>(args)...);
  if (lua_pcall(L, sizeof...(ARGS), 1, 0) == 0)
  {
    if (lua_isnil(L, -1))
      throw std::logic_error("Lua function '" + fn_name + "' did not return anything.");
    ret_val = LuaArgAsType<TRET>(L, -1);
    lua_pop(L, 1);
    return ret_val;
  }
  else
    throw std::logic_error("Call to lua function '" + fn_name + "' failed. ");
}

template <typename T, typename A, typename... ARGS>
inline std::vector<T, A>
LuaCallStdVectorImpl(lua_State* L, const std::string& fn_name, ARGS... args)
{
  std::vector<T, A> ret_val;
  lua_getglobal(L, fn_name.c_str());
  if (not lua_isfunction(L, -1))
    throw std::logic_error("Lua function '" + fn_name + "' does not exist.");

  LuaPushArgs(L, std::forward<ARGS>(args)...);
  if (lua_pcall(L, sizeof...(ARGS), 1, 0) == 0)
  {
    if (lua_isnil(L, -1))
      throw std::logic_error("Lua function '" + fn_name + "' did not return anything.");

    const size_t num = lua_rawlen(L, -1);
    ret_val.reserve(num);
    for (size_t i = 0; i < num; ++i)
    {
      LuaPush(L, i + 1);
      lua_gettable(L, -2);
      ret_val.push_back(LuaArgAsType<T>(L, -1));
      lua_pop(L, 1);
    }
    lua_pop(L, 1);
    return ret_val;
  }
  else
    throw std::logic_error("Call to lua function '" + fn_name + "' failed. ");
}

} // namespace

template <typename TRET, typename... ARGS>
inline TRET
LuaCall(lua_State* L, const std::string& fn_name, ARGS... args)
{
  if constexpr (is_std_vector<TRET>::value)
  {
    using T = typename TRET::value_type;
    return LuaCallStdVectorImpl<T, std::allocator<T>>(L, fn_name, args...);
  }
  else
    return LuaCallImpl<TRET>(L, fn_name, args...);
}

//

template <typename T>
inline void
LuaArgCheckType(const std::string& fn_name, lua_State* L, int index, T val)
{
  throw std::invalid_argument(fn_name + ": Unsupported type: " + std::string(typeid(T).name()));
}

template <>
inline void
LuaArgCheckType(const std::string& fn_name, lua_State* L, int index, char* val)
{
  if (not lua_isstring(L, index))
    throw std::invalid_argument(fn_name + ": Expecting string value for argument " +
                                std::to_string(index));
}

template <>
inline void
LuaArgCheckType(const std::string& fn_name, lua_State* L, int index, const char* val)
{
  if (not lua_isstring(L, index))
    throw std::invalid_argument(fn_name + ": Expecting string value for argument " +
                                std::to_string(index));
}

template <>
inline void
LuaArgCheckType(const std::string& fn_name, lua_State* L, int index, char val)
{
  if (not lua_isinteger(L, index))
    throw std::invalid_argument(fn_name + ": Expecting integer value for argument " +
                                std::to_string(index));
}

template <>
inline void
LuaArgCheckType(const std::string& fn_name, lua_State* L, int index, int val)
{
  if (not lua_isinteger(L, index))
    throw std::invalid_argument(fn_name + ": Expecting integer value for argument " +
                                std::to_string(index));
}

template <>
inline void
LuaArgCheckType(const std::string& fn_name, lua_State* L, int index, long val)
{
  if (not lua_isinteger(L, index))
    throw std::invalid_argument(fn_name + ": Expecting integer value for argument " +
                                std::to_string(index));
}

template <>
inline void
LuaArgCheckType(const std::string& fn_name, lua_State* L, int index, long long val)
{
  if (not lua_isinteger(L, index))
    throw std::invalid_argument(fn_name + ": Expecting integer value for argument " +
                                std::to_string(index));
}

template <>
inline void
LuaArgCheckType(const std::string& fn_name, lua_State* L, int index, unsigned char val)
{
  if (not lua_isinteger(L, index))
    throw std::invalid_argument(fn_name + ": Expecting integer value for argument " +
                                std::to_string(index));
}

template <>
inline void
LuaArgCheckType(const std::string& fn_name, lua_State* L, int index, unsigned int val)
{
  if (not lua_isinteger(L, index))
    throw std::invalid_argument(fn_name + ": Expecting integer value for argument " +
                                std::to_string(index));
}

template <>
inline void
LuaArgCheckType(const std::string& fn_name, lua_State* L, int index, unsigned long val)
{
  if (not lua_isinteger(L, index))
    throw std::invalid_argument(fn_name + ": Expecting integer value for argument " +
                                std::to_string(index));
}

template <>
inline void
LuaArgCheckType(const std::string& fn_name, lua_State* L, int index, unsigned long long val)
{
  if (not lua_isinteger(L, index))
    throw std::invalid_argument(fn_name + ": Expecting integer value for argument " +
                                std::to_string(index));
}

template <>
inline void
LuaArgCheckType(const std::string& fn_name, lua_State* L, int index, float val)
{
  if (not lua_isnumber(L, index))
    throw std::invalid_argument(fn_name + ": Expecting number value for argument " +
                                std::to_string(index));
}

template <>
inline void
LuaArgCheckType(const std::string& fn_name, lua_State* L, int index, double val)
{
  if (not lua_isnumber(L, index))
    throw std::invalid_argument(fn_name + ": Expecting number value for argument " +
                                std::to_string(index));
}

template <>
inline void
LuaArgCheckType(const std::string& fn_name, lua_State* L, int index, bool val)
{
  if (not lua_isboolean(L, index))
    throw std::invalid_argument(fn_name + ": Expecting boolean value for argument " +
                                std::to_string(index));
}

template <>
inline void
LuaArgCheckType(const std::string& fn_name, lua_State* L, int index, std::string val)
{
  if (not lua_isstring(L, index))
    throw std::invalid_argument(fn_name + ": Expecting string value for argument " +
                                std::to_string(index));
}

template <>
inline void
LuaArgCheckType(const std::string& fn_name, lua_State* L, int index, opensn::Vector3 val)
{
  if (not lua_istable(L, index))
    throw std::invalid_argument(fn_name + ": Expecting table value for argument " +
                                std::to_string(index));
}

template <>
inline void
LuaArgCheckType(const std::string& fn_name, lua_State* L, int index, opensn::ParameterBlock val)
{
  if (not lua_istable(L, index))
    throw std::invalid_argument(fn_name + ": Expecting table value for argument " +
                                std::to_string(index));
}

template <typename T, typename A>
inline void
LuaArgCheckStdVectorType(const std::string& fn_name, lua_State* L, int index, std::vector<T, A> val)
{
  if (not lua_istable(L, index))
    throw std::invalid_argument(fn_name + ": Expecting table value for argument " +
                                std::to_string(index));
}

template <class KEY, class VAL, class C, class A>
inline void
LuaArgCheckStdMapType(const std::string& fn_name,
                      lua_State* L,
                      int index,
                      std::map<KEY, VAL, C, A> val)
{
  if (not lua_istable(L, index))
    throw std::invalid_argument(fn_name + ": Expecting table value for argument " +
                                std::to_string(index));
}

/**
 * Check that the lua function argument is of a required type
 *
 * \tparam T C++ type being checked
 * \tparam F Type of the lambda function providing the argument index
 * \param fn_name Lua function name
 * \param L Lua stack
 * \param index Lambda function that provides the index of the argument
 */
template <typename T, typename F>
void
LuaCheckArgWithIndex(const std::string& fn_name, lua_State* L, F&& index)
{
  auto idx = index();
  if constexpr (is_std_vector<T>::value)
  {
    using U = typename T::value_type;
    std::vector<U, std::allocator<U>> val;
    LuaArgCheckStdVectorType<U, std::allocator<U>>(fn_name, L, idx, val);
  }
  else if constexpr (is_std_map<T>::value)
  {
    using KEY = typename T::key_type;
    using VAL = typename T::mapped_type;
    std::map<KEY, VAL> val;
    LuaArgCheckStdMapType<KEY, VAL, std::less<KEY>, std::allocator<std::pair<const KEY, VAL>>>(
      fn_name, L, idx, val);
  }
  else
  {
    T val;
    if constexpr (std::is_arithmetic<T>::value)
      val = 0;
    LuaArgCheckType<T>(fn_name, L, idx, val);
  }
}

/**
 * Check Lua function arguments
 *
 * Use like so:
 * ```
 * LuaCheckArgs<ARG1_TYPE, ARG2_TYPE, ...>(L, "module.LuaFunc");
 * ```
 *
 * \param L Lua stack
 * \param fn_name Name of the Lua function (not the C/C++ function calling this API)
 */
template <typename... ARGS>
inline void
LuaCheckArgs(lua_State* L, const std::string& fn_name)
{
  size_t num_args = lua_gettop(L);
  if (num_args < sizeof...(ARGS))
    throw std::logic_error("Function '" + fn_name + "' expects " + std::to_string(sizeof...(ARGS)) +
                           " arguments.");
  size_t i = 0;
  auto f = [&i]() { return ++i; };
  (LuaCheckArgWithIndex<ARGS>(fn_name, L, f), ...);
}

// API for function return values

/**
 * For functions that return something. This puts all `args` on the Lua stack and returns the
 * number of arguments.
 *
 * Use like this:
 * ```
 * int CFun(lua_State* L) {
 *   return LuaReturn(L, result1, result2, ...);
 * }
 * ```
 *
 * \tparam ARGS Arguments
 * \param L Lua stack
 * \param args Individual arguments that will be put on the Lua stack
 * \return Number of arguments put on Lua stack
 */
template <typename... ARGS>
inline int
LuaReturn(lua_State* L, ARGS&&... args)
{
  LuaPushArgs(L, std::forward<ARGS>(args)...);
  return sizeof...(ARGS);
}

/**
 * For functions that don't have a return value
 *
 * \param L Lua Stack
 * \return Zero since no arguments are put on the Lua stack
 */
inline int
LuaReturn(lua_State* L)
{
  return 0;
}

//

inline int
LuaNumArgs(lua_State* L)
{
  return lua_gettop(L);
}

} // namespace opensnlua
