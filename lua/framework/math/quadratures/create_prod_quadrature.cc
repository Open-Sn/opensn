#include "framework/lua.h"

#include "framework/runtime.h"

#include "framework/math/quadratures/angular_product_quadrature.h"

#include "framework/logging/log.h"

#include <memory>

#include "quadratures_lua.h"
#include "framework/console/console.h"

using namespace opensn;

RegisterLuaFunctionAsIs(chiCreateProductQuadrature);

RegisterLuaConstantAsIs(GAUSS_LEGENDRE, Varying(1));
RegisterLuaConstantAsIs(GAUSS_CHEBYSHEV, Varying(2));
RegisterLuaConstantAsIs(GAUSS_LEGENDRE_LEGENDRE, Varying(3));
RegisterLuaConstantAsIs(GAUSS_LEGENDRE_CHEBYSHEV, Varying(4));
RegisterLuaConstantAsIs(CUSTOM_QUADRATURE, Varying(5));

int
chiCreateProductQuadrature(lua_State* L)
{
  int num_args = lua_gettop(L);
  // Parse argument
  int ident = lua_tonumber(L, 1);
  bool verbose = false;

  if (ident == (int)ProductQuadratureType::GAUSS_LEGENDRE)
  {
    if (num_args < 2) LuaPostArgAmountError("chiCreateProductQuadrature", 2, num_args);

    int Np = lua_tonumber(L, 2);
    if (num_args == 3) verbose = lua_toboolean(L, 3);

    opensn::log.Log() << "Creating Gauss-Legendre Quadrature\n";

    auto new_quad = std::make_shared<AngularQuadratureProdGL>(Np, verbose);

    opensn::Chi::angular_quadrature_stack.push_back(new_quad);
    const size_t index = opensn::Chi::angular_quadrature_stack.size() - 1;
    lua_pushinteger(L, static_cast<lua_Integer>(index));

    if (verbose)
    {
      opensn::log.Log() << "Created Gauss-Legendre Quadrature with " << new_quad->azimu_ang_.size()
                        << " azimuthal angles and " << new_quad->polar_ang_.size()
                        << " polar angles.";
    }

    return 1;
  }
  else if (ident == (int)ProductQuadratureType::GAUSS_LEGENDRE_LEGENDRE)
  {
    if (num_args < 3) LuaPostArgAmountError("chiCreateProductQuadrature", 3, num_args);

    int Na = lua_tonumber(L, 2);
    int Np = lua_tonumber(L, 3);
    if (num_args == 4) verbose = lua_toboolean(L, 4);

    opensn::log.Log() << "Creating Gauss-Legendre-Legendre Quadrature\n";

    auto new_quad = std::make_shared<AngularQuadratureProdGLL>(Na, Np, verbose);

    opensn::Chi::angular_quadrature_stack.push_back(new_quad);
    const size_t index = opensn::Chi::angular_quadrature_stack.size() - 1;
    lua_pushinteger(L, static_cast<lua_Integer>(index));

    if (verbose)
    {
      opensn::log.Log() << "Created Gauss-Legendre-Legendre Quadrature with "
                        << new_quad->azimu_ang_.size() << " azimuthal angles and "
                        << new_quad->polar_ang_.size() << " polar angles.";
    }

    return 1;
  }
  else if (ident == (int)ProductQuadratureType::GAUSS_LEGENDRE_CHEBYSHEV)
  {
    if (num_args < 3) LuaPostArgAmountError("chiCreateProductQuadrature", 3, num_args);

    int Na = lua_tonumber(L, 2);
    int Np = lua_tonumber(L, 3);
    if (num_args == 4) verbose = lua_toboolean(L, 4);

    opensn::log.Log() << "Creating Gauss-Legendre-ChebyShev Quadrature\n";

    auto new_quad = std::make_shared<AngularQuadratureProdGLC>(Na, Np, verbose);

    opensn::Chi::angular_quadrature_stack.push_back(new_quad);
    const size_t index = opensn::Chi::angular_quadrature_stack.size() - 1;
    lua_pushinteger(L, static_cast<lua_Integer>(index));

    if (verbose)
    {
      opensn::log.Log() << "Created Gauss-Legendre-Chebyshev Quadrature with "
                        << new_quad->azimu_ang_.size() << " azimuthal angles and "
                        << new_quad->polar_ang_.size() << " polar angles.";
    }

    return 1;
  }
  else if (ident == (int)ProductQuadratureType::CUSTOM_QUADRATURE)
  {
    if (num_args < 4)
      LuaPostArgAmountError("chiCreateProductQuadrature:CUSTOM_QUADRATURE", 3, num_args);

    if (not lua_istable(L, 2))
    {
      opensn::log.LogAllError()
        << "chiCreateProductQuadrature:CUSTOM_QUADRATURE, second argument must "
        << "be a lua table.";
      opensn::Exit(EXIT_FAILURE);
    }
    if (not lua_istable(L, 3))
    {
      opensn::log.LogAllError()
        << "chiCreateProductQuadrature:CUSTOM_QUADRATURE, third argument must "
        << "be a lua table.";
      opensn::Exit(EXIT_FAILURE);
    }
    if (not lua_istable(L, 4))
    {
      opensn::log.LogAllError()
        << "chiCreateProductQuadrature:CUSTOM_QUADRATURE, fourth argument must "
        << "be a lua table.";
      opensn::Exit(EXIT_FAILURE);
    }
    if (num_args == 5) verbose = lua_toboolean(L, 4);

    size_t Na = lua_rawlen(L, 2);
    size_t Np = lua_rawlen(L, 3);
    size_t Nw = lua_rawlen(L, 4);

    std::vector<double> azimuthal(Na, 0.0);
    std::vector<double> polar(Np, 0.0);
    std::vector<double> weights(Nw, 0.0);

    for (int n = 1; n <= Na; ++n)
    {
      lua_pushnumber(L, n);
      lua_gettable(L, 2);
      azimuthal[n - 1] = lua_tonumber(L, -1);
      lua_pop(L, 1);
    }
    for (int n = 1; n <= Np; ++n)
    {
      lua_pushnumber(L, n);
      lua_gettable(L, 3);
      polar[n - 1] = lua_tonumber(L, -1);
      lua_pop(L, 1);
    }
    for (int n = 1; n <= Nw; ++n)
    {
      lua_pushnumber(L, n);
      lua_gettable(L, 4);
      weights[n - 1] = lua_tonumber(L, -1);
      lua_pop(L, 1);
    }

    opensn::log.Log() << "Creating custom product quadrature Quadrature\n";

    opensn::log.Log() << Na << " " << Np << " " << Nw;

    auto new_quad =
      std::make_shared<AngularQuadratureProdCustom>(azimuthal, polar, weights, verbose);

    opensn::Chi::angular_quadrature_stack.push_back(new_quad);
    const size_t index = opensn::Chi::angular_quadrature_stack.size() - 1;
    lua_pushinteger(L, static_cast<lua_Integer>(index));

    if (verbose)
    {
      opensn::log.Log() << "Created Custom Quadrature with " << new_quad->azimu_ang_.size()
                        << " azimuthal angles and " << new_quad->polar_ang_.size()
                        << " polar angles.";
    }

    return 1;
  }
  else
  {
    opensn::log.LogAllError()
      << "In call to chiCreateProductQuadrature. Unsupported quadrature type"
         " supplied. Given: "
      << ident;
    opensn::Exit(EXIT_FAILURE);
  }
  return 0;
}
