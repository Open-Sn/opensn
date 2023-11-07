#include "framework/lua.h"

#include "framework/runtime.h"

#include "framework/math/quadratures/angular_quadrature_base.h"

#include "framework/logging/log.h"

#include "quadratures_lua.h"
#include "framework/console/console.h"

RegisterLuaFunctionAsIs(chiCreateCustomAngularQuadrature);

int
chiCreateCustomAngularQuadrature(lua_State* L)
{
  const std::string fname = __FUNCTION__;
  const int num_args = lua_gettop(L);

  if (num_args != 3) LuaPostArgAmountError(fname, 3, num_args);

  LuaCheckNilValue(fname, L, 1);
  LuaCheckNilValue(fname, L, 2);
  LuaCheckNilValue(fname, L, 3);

  LuaCheckTableValue(fname, L, 1);
  LuaCheckTableValue(fname, L, 2);
  LuaCheckTableValue(fname, L, 3);

  size_t Na = lua_rawlen(L, 1);
  size_t Np = lua_rawlen(L, 2);
  size_t Nw = lua_rawlen(L, 3);

  if ((Na - Np != 0) or (Na - Nw != 0))
  {
    Chi::log.LogAllError() << fname + ": Tables lengths supplied "
                                      "are not of equal lengths.";
    Chi::Exit(EXIT_FAILURE);
  }

  std::vector<double> azi_angles(Na, 0.0);
  std::vector<double> pol_angles(Na, 0.0);
  std::vector<double> weights(Na, 0.0);

  for (int n = 1; n <= Na; ++n)
  {
    lua_pushnumber(L, n);
    lua_gettable(L, 1);
    azi_angles[n - 1] = lua_tonumber(L, -1);
    lua_pop(L, 1);
  }

  for (int n = 1; n <= Na; ++n)
  {
    lua_pushnumber(L, n);
    lua_gettable(L, 2);
    pol_angles[n - 1] = lua_tonumber(L, -1);
    lua_pop(L, 1);
  }

  for (int n = 1; n <= Na; ++n)
  {
    lua_pushnumber(L, n);
    lua_gettable(L, 3);
    weights[n - 1] = lua_tonumber(L, -1);
    lua_pop(L, 1);
  }

  Chi::log.Log() << "Creating Custom Angular Quadrature\n";

  auto new_quad =
    std::make_shared<chi_math::AngularQuadratureCustom>(azi_angles, pol_angles, weights, false);

  Chi::angular_quadrature_stack.push_back(new_quad);
  size_t index = Chi::angular_quadrature_stack.size() - 1;
  lua_pushinteger(L, static_cast<lua_Integer>(index));

  return 1;
}
