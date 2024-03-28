#include "framework/lua.h"
#include "framework/math/quadratures/angular/legendre_poly/legendrepoly.h"
#include "framework/console/console.h"
#include "legendre_lua.h"

using namespace opensn;

namespace opensnlua
{

RegisterLuaFunctionNamespace(Legendre, aquad, Legendre);
RegisterLuaFunctionNamespace(LegendreDerivative, aquad, LegendreDerivative);
RegisterLuaFunctionNamespace(Ylm, aquad, Ylm);

int
Legendre(lua_State* L)
{
  // Retrieve arguments
  int N = lua_tonumber(L, 1);
  double x = lua_tonumber(L, 2);

  double retval = opensn::Legendre(N, x);

  lua_pushnumber(L, retval);
  return 1;
}

int
LegendreDerivative(lua_State* L)
{
  // Retrieve arguments
  int N = lua_tonumber(L, 1);
  double x = lua_tonumber(L, 2);

  double retval = dLegendredx(N, x);

  lua_pushnumber(L, retval);
  return 1;
}

int
Ylm(lua_State* L)
{
  int num_args = lua_gettop(L);
  if (num_args != 4)
    LuaPostArgAmountError("Ylm", 4, num_args);

  int ell = lua_tonumber(L, 1);
  int m = lua_tonumber(L, 2);
  double theta = lua_tonumber(L, 3);
  double varphi = lua_tonumber(L, 4);

  double retval = opensn::Ylm(ell, m, varphi, theta);

  lua_pushnumber(L, retval);
  return 1;
}

} // namespace opensnlua
