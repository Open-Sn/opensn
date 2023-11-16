#include "framework/lua.h"

#include "framework/runtime.h"

#include "framework/math/quadratures/angular_product_quadrature.h"

#include "framework/logging/log.h"

#include "quadratures_lua.h"
#include "framework/console/console.h"

using namespace opensn;

RegisterLuaFunctionAsIs(chiGetProductQuadrature);

int
chiGetProductQuadrature(lua_State* L)
{
  int num_args = lua_gettop(L);
  if (num_args != 1) LuaPostArgAmountError("chiGetProductQuadrature", 1, num_args);

  int handle = lua_tonumber(L, 1);

  std::shared_ptr<ProductQuadrature> quad;
  try
  {
    auto ang_quad = opensn::Chi::angular_quadrature_stack.at(handle);
    if (ang_quad->type_ == AngularQuadratureType::ProductQuadrature)
      quad = std::static_pointer_cast<ProductQuadrature>(ang_quad);
    else
    {
      opensn::log.LogAllError() << "chiGetProductQuadrature: Provided quadrature handle points to "
                                   "a quadrature that is not a product quadrature.";
      opensn::Exit(EXIT_FAILURE);
    }
  }
  catch (const std::out_of_range& o)
  {
    opensn::log.LogAllError() << "chiGetProductQuadrature: Invalid quadrature handle.";
    opensn::Exit(EXIT_FAILURE);
  }

  lua_newtable(L);
  for (size_t n = 0; n < quad->weights_.size(); ++n)
  {
    lua_pushnumber(L, n + 1);
    lua_newtable(L);

    lua_pushstring(L, "weight");
    lua_pushnumber(L, quad->weights_[n]);
    lua_settable(L, -3);

    lua_pushstring(L, "polar");
    lua_pushnumber(L, quad->abscissae_[n].theta);
    lua_settable(L, -3);

    lua_pushstring(L, "azimuthal");
    lua_pushnumber(L, quad->abscissae_[n].phi);
    lua_settable(L, -3);

    lua_settable(L, -3);
  }

  return 1;
}
