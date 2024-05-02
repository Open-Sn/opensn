// SPDX-FileCopyrightText: 2024 The OpenSn Authors <https://open-sn.github.io/opensn/>
// SPDX-License-Identifier: MIT

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

RegisterLuaConstant(GAUSS_LEGENDRE, Varying(1));
RegisterLuaConstant(GAUSS_CHEBYSHEV, Varying(2));
RegisterLuaConstant(GAUSS_LEGENDRE_LEGENDRE, Varying(3));
RegisterLuaConstant(GAUSS_LEGENDRE_CHEBYSHEV, Varying(4));
RegisterLuaConstant(CUSTOM_QUADRATURE, Varying(5));

int
CreateProductQuadrature(lua_State* L)
{
  const std::string fname = "aquad.CreateProductQuadrature";

  auto ident = LuaArg<int>(L, 1);
  if (ident == (int)ProductQuadratureType::GAUSS_LEGENDRE)
  {
    LuaCheckArgs<int, int>(L, fname);

    auto Np = LuaArg<int>(L, 2);
    auto verbose = LuaArgOptional<bool>(L, 3, false);

    opensn::log.Log() << "Creating Gauss-Legendre Quadrature\n";

    auto new_quad = std::make_shared<AngularQuadratureProdGL>(Np, verbose);

    opensn::angular_quadrature_stack.push_back(new_quad);
    const size_t index = opensn::angular_quadrature_stack.size() - 1;

    if (verbose)
    {
      opensn::log.Log() << "Created Gauss-Legendre Quadrature with " << new_quad->azimu_ang_.size()
                        << " azimuthal angles and " << new_quad->polar_ang_.size()
                        << " polar angles.";
    }

    return LuaReturn(L, index);
  }
  else if (ident == (int)ProductQuadratureType::GAUSS_LEGENDRE_LEGENDRE)
  {
    LuaCheckArgs<int, int, int>(L, fname);

    auto Na = LuaArg<int>(L, 2);
    auto Np = LuaArg<int>(L, 3);
    auto verbose = LuaArgOptional<bool>(L, 4, false);

    opensn::log.Log() << "Creating Gauss-Legendre-Legendre Quadrature\n";

    auto new_quad = std::make_shared<AngularQuadratureProdGLL>(Na, Np, verbose);

    opensn::angular_quadrature_stack.push_back(new_quad);
    const size_t index = opensn::angular_quadrature_stack.size() - 1;

    if (verbose)
    {
      opensn::log.Log() << "Created Gauss-Legendre-Legendre Quadrature with "
                        << new_quad->azimu_ang_.size() << " azimuthal angles and "
                        << new_quad->polar_ang_.size() << " polar angles.";
    }

    return LuaReturn(L, index);
  }
  else if (ident == (int)ProductQuadratureType::GAUSS_LEGENDRE_CHEBYSHEV)
  {
    LuaCheckArgs<int, int, int>(L, fname);

    auto Na = LuaArg<int>(L, 2);
    auto Np = LuaArg<int>(L, 3);
    auto verbose = LuaArgOptional<bool>(L, 4, false);

    opensn::log.Log() << "Creating Gauss-Legendre-ChebyShev Quadrature\n";

    auto new_quad = std::make_shared<AngularQuadratureProdGLC>(Na, Np, verbose);

    opensn::angular_quadrature_stack.push_back(new_quad);
    const size_t index = opensn::angular_quadrature_stack.size() - 1;

    if (verbose)
    {
      opensn::log.Log() << "Created Gauss-Legendre-Chebyshev Quadrature with "
                        << new_quad->azimu_ang_.size() << " azimuthal angles and "
                        << new_quad->polar_ang_.size() << " polar angles.";
    }

    return LuaReturn(L, index);
  }
  else if (ident == (int)ProductQuadratureType::CUSTOM_QUADRATURE)
  {
    LuaCheckArgs<int, std::vector<double>, std::vector<double>, std::vector<double>>(L, fname);

    auto azimuthal = LuaArg<std::vector<double>>(L, 2);
    auto polar = LuaArg<std::vector<double>>(L, 3);
    auto weights = LuaArg<std::vector<double>>(L, 4);
    auto verbose = LuaArgOptional<bool>(L, 5, false);

    opensn::log.Log() << "Creating custom product quadrature Quadrature\n";

    opensn::log.Log() << azimuthal.size() << " " << polar.size() << " " << weights.size();

    auto new_quad =
      std::make_shared<AngularQuadratureProdCustom>(azimuthal, polar, weights, verbose);

    opensn::angular_quadrature_stack.push_back(new_quad);
    const size_t index = opensn::angular_quadrature_stack.size() - 1;

    if (verbose)
    {
      opensn::log.Log() << "Created Custom Quadrature with " << new_quad->azimu_ang_.size()
                        << " azimuthal angles and " << new_quad->polar_ang_.size()
                        << " polar angles.";
    }

    return LuaReturn(L, index);
  }
  else
  {
    opensn::log.LogAllError() << "In call to CreateProductQuadrature. Unsupported quadrature type"
                                 " supplied. Given: "
                              << ident;
    opensn::Exit(EXIT_FAILURE);
  }
  return LuaReturn(L);
}

} // namespace opensnlua
