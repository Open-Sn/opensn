// SPDX-FileCopyrightText: 2024 The OpenSn Authors <https://open-sn.github.io/opensn/>
// SPDX-License-Identifier: MIT

#include "lua/framework/lua.h"
#include "framework/math/quadratures/angular/legendre_poly/legendrepoly.h"
#include "lua/framework/console/console.h"
#include "lua/framework/math/quadratures/legendre_poly/legendre.h"

using namespace opensn;

namespace opensnlua
{

RegisterLuaFunctionInNamespace(Legendre, aquad, Legendre);
RegisterLuaFunctionInNamespace(LegendreDerivative, aquad, LegendreDerivative);
RegisterLuaFunctionInNamespace(Ylm, aquad, Ylm);

int
Legendre(lua_State* L)
{
  LuaCheckArgs<int, double>(L, "aquad.Legendre");
  auto N = LuaArg<int>(L, 1);
  auto x = LuaArg<double>(L, 2);
  double retval = opensn::Legendre(N, x);
  return LuaReturn(L, retval);
}

int
LegendreDerivative(lua_State* L)
{
  LuaCheckArgs<int, double>(L, "aquad.LegendreDerivative");
  auto N = LuaArg<int>(L, 1);
  auto x = LuaArg<double>(L, 2);
  double retval = dLegendredx(N, x);
  return LuaReturn(L, retval);
}

int
Ylm(lua_State* L)
{
  LuaCheckArgs<int, int, double, double>(L, "aquad.Ylm");
  auto ell = LuaArg<int>(L, 1);
  auto m = LuaArg<int>(L, 2);
  auto theta = LuaArg<double>(L, 3);
  auto varphi = LuaArg<double>(L, 4);
  double retval = opensn::Ylm(ell, m, varphi, theta);
  return LuaReturn(L, retval);
}

} // namespace opensnlua
