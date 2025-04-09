// SPDX-FileCopyrightText: 2024 The OpenSn Authors <https://open-sn.github.io/opensn/>
// SPDX-License-Identifier: MIT

#pragma once

#include "framework/parameters/parameter_block.h"
#include "lua/lib/parse_table.h"
#include "framework/parameters/input_parameters.h"
#include "framework/math/vector3.h"
#include "framework/math/quadratures/angular/angular_quadrature.h"
#include "LuaBridge/LuaBridge.h"
#include <stdexcept>
extern "C"
{
#include "lua.h"
#include "lualib.h"
#include "lauxlib.h"
}
#include <string>

namespace luabridge
{

template <>
struct Stack<opensn::InputParameters>
{
  inline static opensn::InputParameters get(lua_State* L, int idx)
  {
    opensn::InputParameters result;
    opensnlua::ParseTableKeys(L, idx, result);
    return result;
  }

  inline static Result push(lua_State* L, const opensn::InputParameters& v)
  {
    lua_newtable(L);
    // TODO: put params on the stack
    return {};
  }

  inline static bool isInstance(lua_State* L, int index)
  {
    return lua_type(L, index) == LUA_TTABLE;
  }
};

template <>
struct Stack<opensn::ParameterBlock>
{
  inline static TypeResult<opensn::ParameterBlock> get(lua_State* L, int idx)
  {
    opensn::ParameterBlock result;
    opensnlua::ParseTableKeys(L, idx, result);
    return result;
  }

  inline static Result push(lua_State* L, const opensn::ParameterBlock& v)
  {
    return opensnlua::PushParameterBlock(L, v, 0);
  }

  inline static bool isInstance(lua_State* L, int index)
  {
    return lua_type(L, index) == LUA_TTABLE;
  }
};

template <>
struct Stack<opensn::Vector3>
{
  inline static TypeResult<opensn::Vector3> get(lua_State* L, int idx)
  {
    opensn::Vector3 vec;
    auto res = Stack<std::string>::push(L, "x");
    lua_gettable(L, idx);
    vec.x = Stack<double>::get(L, -1).value();
    lua_pop(L, 1);

    res = Stack<std::string>::push(L, "y");
    lua_gettable(L, idx);
    if (not lua_isnil(L, -1))
      vec.y = Stack<double>::get(L, -1).value();
    lua_pop(L, 1);

    res = Stack<std::string>::push(L, "z");
    lua_gettable(L, idx);
    if (not lua_isnil(L, -1))
      vec.z = Stack<double>::get(L, -1).value();
    lua_pop(L, 1);

    return vec;
  }

  inline static Result push(lua_State* L, const opensn::Vector3& v)
  {
    lua_newtable(L);
    auto res = Stack<std::string>::push(L, "x");
    res = Stack<double>::push(L, v.x);
    lua_settable(L, -3);
    res = Stack<std::string>::push(L, "y");
    res = Stack<double>::push(L, v.y);
    lua_settable(L, -3);
    res = Stack<std::string>::push(L, "z");
    res = Stack<double>::push(L, v.z);
    lua_settable(L, -3);
    return {};
  }

  inline static bool isInstance(lua_State* L, int index)
  {
    return lua_type(L, index) == LUA_TTABLE;
  }
};

template <>
struct Stack<opensn::QuadraturePointPhiTheta>
{
  inline static TypeResult<opensn::QuadraturePointPhiTheta> get(lua_State* L, int idx)
  {
    auto res = Stack<std::string>::push(L, "phi");
    lua_gettable(L, idx);
    auto phi = Stack<double>::get(L, -1).value();
    lua_pop(L, 1);

    res = Stack<std::string>::push(L, "theta");
    lua_gettable(L, idx);
    auto theta = Stack<double>::get(L, -1).value();
    lua_pop(L, 1);

    return opensn::QuadraturePointPhiTheta(phi, theta);
  }

  inline static Result push(lua_State* L, const opensn::QuadraturePointPhiTheta& pt)
  {
    lua_newtable(L);
    auto res = Stack<std::string>::push(L, "phi");
    res = Stack<double>::push(L, pt.phi);
    lua_settable(L, -3);
    res = Stack<std::string>::push(L, "theta");
    res = Stack<double>::push(L, pt.theta);
    lua_settable(L, -3);
    return {};
  }

  inline static bool isInstance(lua_State* L, int index)
  {
    return lua_type(L, index) == LUA_TTABLE;
  }
};

} // namespace luabridge
