#include "framework/chi_lua.h"

#include "framework/chi_runtime.h"
#include "framework/logging/chi_log.h"

#include "framework/ChiObjectFactory.h"

#include "quadratures_lua.h"
#include "framework/console/chi_console.h"

RegisterLuaFunctionAsIs(chiCreateLineQuadrature);

int
chiCreateLineQuadrature(lua_State* L)
{
  const std::string fname = __FUNCTION__;
  const int num_args = lua_gettop(L);

  if (not((num_args == 2) or (num_args == 3))) LuaPostArgAmountError(fname, 2, num_args);

  LuaCheckNilValue(fname, L, 1);
  LuaCheckNilValue(fname, L, 2);

  // Parse argument
  int ident = lua_tonumber(L, 1);
  int N = lua_tonumber(L, 2);
  bool verbose = false;
  if (num_args == 3) verbose = lua_toboolean(L, 3);

  chi::ParameterBlock params;
  params.AddParameter("verbose", verbose);
  params.AddParameter("N", N);

  auto& obj_factory = ChiObjectFactory::GetInstance();

  if (ident == 1) // GAUSS_LEGENDRE
  {
    Chi::log.Log() << "Creating Gauss-Legendre Quadrature\n";

    const size_t handle =
      obj_factory.MakeRegisteredObjectOfType("chi_math::QuadratureGaussLegendre", params);

    lua_pushinteger(L, static_cast<lua_Integer>(handle));
    return 1;
  }
  else if (ident == 2) // GAUSS_CHEBYSHEV
  {
    Chi::log.Log() << "Creating Gauss-Chebyshev Quadrature\n";

    const size_t handle =
      obj_factory.MakeRegisteredObjectOfType("chi_math::QuadratureGaussChebyshev", params);

    lua_pushinteger(L, static_cast<lua_Integer>(handle));
    return 1;
  }
  return 0;
}
