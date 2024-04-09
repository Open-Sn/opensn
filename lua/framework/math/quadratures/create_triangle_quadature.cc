// SPDX-FileCopyrightText: 2024 The OpenSn Authors <https://open-sn.github.io/opensn/>
// SPDX-License-Identifier: MIT

#include "framework/lua.h"
#include "quadratures_lua.h"
#include "framework/math/quadratures/angular/triangle_quadrature.h"
#include "framework/logging/log.h"
#include "framework/console/console.h"
#include "framework/runtime.h"
#include <memory>

using namespace opensn;

namespace opensnlua
{

RegisterLuaFunctionNamespace(CreateTriangleQuadrature, aquad, CreateTriangleQuadrature);

int
CreateTriangleQuadrature(lua_State* L)
{
  int num_args = lua_gettop(L);
  // Parse argument
  int method = lua_tonumber(L, 1);
  int order = lua_tonumber(L, 2);
  int moments = 0;
  if (num_args == 3)
    moments = lua_tonumber(L, 3);

  if (num_args < 2)
    LuaPostArgAmountError("CreateTriangleQuadrature", 2, num_args);

  opensn::log.Log() << "Creating Triangle Quadrature\n";

  auto new_quad = std::make_shared<TriangleQuadrature>(method, order, moments);
  opensn::angular_quadrature_stack.push_back(new_quad);
  const size_t index = opensn::angular_quadrature_stack.size() - 1;
  lua_pushinteger(L, static_cast<lua_Integer>(index));

  opensn::log.Log() << "Created Triangle Quadrature with method = " << method << " and order "
                    << order << std::endl;

  return 1;
}

} // namespace opensnlua
