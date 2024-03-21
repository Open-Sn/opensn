#include "framework/lua.h"
#include "quadratures_lua.h"
#include "framework/math/quadratures/angular/product_quadrature.h"
#include "framework/logging/log.h"
#include "framework/console/console.h"
#include "framework/runtime.h"
#include <memory>

using namespace opensn;

namespace opensnlua
{

RegisterLuaFunctionNamespace(CreateProductQuadrature, aquad, CreateProductQuadrature);

RegisterLuaConstantAsIs(GAUSS_LEGENDRE, Varying(1));
RegisterLuaConstantAsIs(GAUSS_CHEBYSHEV, Varying(2));
RegisterLuaConstantAsIs(GAUSS_LEGENDRE_LEGENDRE, Varying(3));
RegisterLuaConstantAsIs(GAUSS_LEGENDRE_CHEBYSHEV, Varying(4));
RegisterLuaConstantAsIs(CUSTOM_QUADRATURE, Varying(5));

int
CreateProductQuadrature(lua_State* L)
{
  int num_args = lua_gettop(L);
  // Parse argument
  auto ident = LuaArg<int>(L, 1);

  if (ident == (int)ProductQuadratureType::GAUSS_LEGENDRE)
  {
    if (num_args < 2)
      LuaPostArgAmountError("CreateProductQuadrature", 2, num_args);

    auto Np = LuaArg<int>(L, 2);
    auto verbose = LuaArgOptional<bool>(L, 3, false);

    opensn::log.Log() << "Creating Gauss-Legendre Quadrature\n";

    auto new_quad = std::make_shared<AngularQuadratureProdGL>(Np, verbose);

    opensn::angular_quadrature_stack.push_back(new_quad);
    const size_t index = opensn::angular_quadrature_stack.size() - 1;
    LuaPush(L, index);

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
    if (num_args < 3)
      LuaPostArgAmountError("CreateProductQuadrature", 3, num_args);

    auto Na = LuaArg<int>(L, 2);
    auto Np = LuaArg<int>(L, 3);
    auto verbose = LuaArgOptional<bool>(L, 4, false);

    opensn::log.Log() << "Creating Gauss-Legendre-Legendre Quadrature\n";

    auto new_quad = std::make_shared<AngularQuadratureProdGLL>(Na, Np, verbose);

    opensn::angular_quadrature_stack.push_back(new_quad);
    const size_t index = opensn::angular_quadrature_stack.size() - 1;
    LuaPush(L, index);

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
    if (num_args < 3)
      LuaPostArgAmountError("CreateProductQuadrature", 3, num_args);

    auto Na = LuaArg<int>(L, 2);
    auto Np = LuaArg<int>(L, 3);
    auto verbose = LuaArgOptional<bool>(L, 4, false);

    opensn::log.Log() << "Creating Gauss-Legendre-ChebyShev Quadrature\n";

    auto new_quad = std::make_shared<AngularQuadratureProdGLC>(Na, Np, verbose);

    opensn::angular_quadrature_stack.push_back(new_quad);
    const size_t index = opensn::angular_quadrature_stack.size() - 1;
    LuaPush(L, index);

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
      LuaPostArgAmountError("CreateProductQuadrature:CUSTOM_QUADRATURE", 3, num_args);

    std::vector<double> azimuthal = LuaArgVector<double>(L, 2);
    std::vector<double> polar = LuaArgVector<double>(L, 3);
    std::vector<double> weights = LuaArgVector<double>(L, 4);
    auto verbose = LuaArgOptional<bool>(L, 5, false);

    opensn::log.Log() << "Creating custom product quadrature Quadrature\n";

    opensn::log.Log() << azimuthal.size() << " " << polar.size() << " " << weights.size();

    auto new_quad =
      std::make_shared<AngularQuadratureProdCustom>(azimuthal, polar, weights, verbose);

    opensn::angular_quadrature_stack.push_back(new_quad);
    const size_t index = opensn::angular_quadrature_stack.size() - 1;
    LuaPush(L, index);

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
    opensn::log.LogAllError() << "In call to CreateProductQuadrature. Unsupported quadrature type"
                                 " supplied. Given: "
                              << ident;
    opensn::Exit(EXIT_FAILURE);
  }
  return 0;
}

} // namespace opensnlua
