#include "framework/lua.h"
#include "quadratures_lua.h"
#include "framework/math/quadratures/angular/product_quadrature.h"
#include "framework/logging/log.h"
#include "framework/console/console.h"
#include "framework/runtime.h"

using namespace opensn;

namespace opensnlua
{

RegisterLuaFunctionNamespace(GetProductQuadrature, aquad, GetProductQuadrature);

int
GetProductQuadrature(lua_State* L)
{
  const std::string fname = "aquad.GetProductQuadrature";
  LuaCheckArgs<int>(L, fname);

  auto handle = LuaArg<size_t>(L, 1);

  std::shared_ptr<ProductQuadrature> quad;
  try
  {
    auto ang_quad = opensn::angular_quadrature_stack.at(handle);
    if (ang_quad->type_ == AngularQuadratureType::ProductQuadrature)
      quad = std::static_pointer_cast<ProductQuadrature>(ang_quad);
    else
    {
      opensn::log.LogAllError() << fname + ": Provided quadrature handle points to a quadrature "
                                           "that is not a product quadrature.";
      opensn::Exit(EXIT_FAILURE);
    }
  }
  catch (const std::out_of_range& o)
  {
    opensn::log.LogAllError() << fname + ": Invalid quadrature handle.";
    opensn::Exit(EXIT_FAILURE);
  }

  lua_newtable(L);
  for (size_t n = 0; n < quad->weights_.size(); ++n)
  {
    LuaPush(L, n + 1);
    lua_newtable(L);
    LuaPushTableKey(L, "weight", quad->weights_[n]);
    LuaPushTableKey(L, "polar", quad->abscissae_[n].theta);
    LuaPushTableKey(L, "azimuthal", quad->abscissae_[n].phi);
    lua_settable(L, -3);
  }

  return 1;
}

} // namespace opensnlua
