#include "framework/lua.h"

#include "framework/runtime.h"

#include "framework/math/quadratures/quadrature_gausschebyshev.h"
#include "framework/math/quadratures/quadrature_gausslegendre.h"
#include "framework/math/quadratures/cylindrical_angular_quadrature.h"
#include "framework/math/quadratures/spherical_angular_quadrature.h"

#include "framework/logging/log.h"

#include "quadratures_lua.h"
#include "framework/console/console.h"

using namespace opensn;

RegisterLuaFunctionAsIs(chiCreateCylindricalProductQuadrature);
RegisterLuaFunctionAsIs(chiCreateSphericalProductQuadrature);

int
chiCreateCylindricalProductQuadrature(lua_State* L)
{
  const std::string fname = __FUNCTION__;
  const int num_args = lua_gettop(L);
  if (num_args < 3) LuaPostArgAmountError(fname, 3, num_args);

  const int ident = lua_tonumber(L, 1);
  const int Np = lua_tonumber(L, 2);

  std::vector<int> vNa;
  if (lua_isnumber(L, 3))
  {
    const int Na = lua_tonumber(L, 3);
    vNa.resize(Np, Na);
  }
  else if (lua_istable(L, 3))
  {
    const size_t lNa = lua_rawlen(L, 3);
    if (lNa != Np)
    {
      opensn::Chi::log.LogAllError()
        << "chiCreateCylindricalProductQuadrature : third argument, "
        << ", if a lua table, must be of length equal to second argument.";
      opensn::Chi::Exit(EXIT_FAILURE);
    }
    vNa.resize(Np, 0);
    for (int n = 1; n <= lNa; ++n)
    {
      lua_pushnumber(L, n);
      lua_gettable(L, 3);
      vNa[n - 1] = lua_tonumber(L, -1);
      lua_pop(L, 1);
    }
  }
  else
  {
    opensn::Chi::log.LogAllError() << "chiCreateCylindricalProductQuadrature : third argument "
                                   << "must be a number or a lua table.";
    opensn::Chi::Exit(EXIT_FAILURE);
  }

  bool verbose = false;
  if (num_args == 4) verbose = lua_toboolean(L, 4);

  const auto prod_quad_type = static_cast<ProductQuadratureType>(ident);
  switch (prod_quad_type)
  {
    case ProductQuadratureType::GAUSS_LEGENDRE_CHEBYSHEV:
    {
      opensn::Chi::log.Log() << "chiCreateCylindricalProductQuadrature : "
                             << "Creating Gauss-Legendre-Legendre Quadrature\n";

      const auto quad_pol = QuadratureGaussLegendre(Np, verbose);
      std::vector<Quadrature> quad_azi;
      for (const auto& Na : vNa)
        quad_azi.emplace_back(QuadratureGaussChebyshev(Na, verbose));
      const auto new_quad =
        std::make_shared<CylindricalAngularQuadrature>(quad_pol, quad_azi, verbose);

      opensn::Chi::angular_quadrature_stack.push_back(new_quad);
      const size_t index = opensn::Chi::angular_quadrature_stack.size() - 1;
      lua_pushinteger(L, static_cast<lua_Integer>(index));

      return 1;
    }
    case ProductQuadratureType::GAUSS_LEGENDRE_LEGENDRE:
    {
      opensn::Chi::log.Log() << "chiCreateCylindricalProductQuadrature : "
                             << "Creating Gauss-Legendre-Legendre Quadrature\n";

      const auto quad_pol = QuadratureGaussLegendre(Np, verbose);
      std::vector<Quadrature> quad_azi;
      for (const auto& Na : vNa)
        quad_azi.emplace_back(QuadratureGaussLegendre(Na, verbose));
      const auto new_quad =
        std::make_shared<CylindricalAngularQuadrature>(quad_pol, quad_azi, verbose);

      opensn::Chi::angular_quadrature_stack.push_back(new_quad);
      const size_t index = opensn::Chi::angular_quadrature_stack.size() - 1;
      lua_pushinteger(L, static_cast<lua_Integer>(index));

      return 1;
    }
    default:
    {
      opensn::Chi::log.LogAllError() << "chiCreateCylindricalProductQuadrature : "
                                     << "Unsupported quadrature type supplied, type=" << ident;
      opensn::Chi::Exit(EXIT_FAILURE);
    }
  }

  return 0;
}

/** Creates a curvilinear product quadrature suitable for spherical geometries.

 \param QuadratureType int Quadrature identifier.
 \param values varying Varying options based on the quadrature type.

 ##_

 ###QuadratureType:\n
 GAUSS_CHEBYSHEV\n
   Gauss-Chebyshev quadrature for the polar angle and no quadrature
   for the azimuthal angle.
   Arguments for this quadrature type, in order:
   - Np : (int) number of polar angles
   - verbose : (bool) verbosity flag (optional).

 ###QuadratureType:\n
 GAUSS_LEGENDRE\n
   Gauss-Legendre quadrature for the polar angle and no quadrature
   for the azimuthal angle.
   Arguments for this quadrature type, in order:
   - Np : (int) number of polar angles
   - verbose : (bool) verbosity flag (optional).


 \return Returns a unique handle to the created product quadrature rule

 \ingroup LuaQuadrature
 */
int
chiCreateSphericalProductQuadrature(lua_State* L)
{
  const int num_args = lua_gettop(L);
  if (num_args < 2) LuaPostArgAmountError("chiCreateSphericalProductQuadrature", 2, num_args);

  const int ident = lua_tonumber(L, 1);
  const int Np = lua_tonumber(L, 2);
  bool verbose = false;
  if (num_args == 3) verbose = lua_toboolean(L, 3);

  const auto prod_quad_type = static_cast<ProductQuadratureType>(ident);
  switch (prod_quad_type)
  {
    case ProductQuadratureType::GAUSS_CHEBYSHEV:
    {
      opensn::Chi::log.Log() << "chiCreateSphericalProductQuadrature : "
                             << "Creating Gauss-Chebyshev Quadrature\n";

      const auto quad_pol = QuadratureGaussChebyshev(Np, verbose);
      const auto new_quad = std::make_shared<SphericalAngularQuadrature>(quad_pol, verbose);

      opensn::Chi::angular_quadrature_stack.push_back(new_quad);
      const size_t index = opensn::Chi::angular_quadrature_stack.size() - 1;
      lua_pushnumber(L, static_cast<lua_Number>(index));

      return 1;
    }
    case ProductQuadratureType::GAUSS_LEGENDRE:
    {
      opensn::Chi::log.Log() << "chiCreateSphericalProductQuadrature : "
                             << "Creating Gauss-Legendre Quadrature\n";

      const auto quad_pol = QuadratureGaussLegendre(Np, verbose);
      const auto new_quad = std::make_shared<SphericalAngularQuadrature>(quad_pol, verbose);

      opensn::Chi::angular_quadrature_stack.push_back(new_quad);
      const size_t index = opensn::Chi::angular_quadrature_stack.size() - 1;
      lua_pushnumber(L, static_cast<lua_Number>(index));

      return 1;
    }
    default:
    {
      opensn::Chi::log.LogAllError() << "chiCreateSphericalProductQuadrature : "
                                     << "Unsupported quadrature type supplied, type=" << ident;
      opensn::Chi::Exit(EXIT_FAILURE);
    }
  }

  return 0;
}
