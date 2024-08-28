// SPDX-FileCopyrightText: 2024 The OpenSn Authors <https://open-sn.github.io/opensn/>
// SPDX-License-Identifier: MIT

#include "lua/framework/lua.h"
#include "lua/framework/math/quadratures/quadratures.h"
#include "lua/framework/console/console.h"
#include "framework/math/quadratures/angular/product_quadrature.h"
#include "framework/logging/log.h"
#include "framework/runtime.h"

using namespace opensn;

namespace opensnlua
{

RegisterLuaFunctionInNamespace(GetProductQuadrature, aquad, GetProductQuadrature);

struct LuaQuadData
{
  double weight;
  double phi;
  double theta;
};

void
LuaPush(lua_State* L, const LuaQuadData& data)
{
  lua_newtable(L);
  LuaPushTableKey(L, "weight", data.weight);
  LuaPushTableKey(L, "polar", data.theta);
  LuaPushTableKey(L, "azimuthal", data.phi);
}

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

  std::vector<LuaQuadData> quad_data;
  quad_data.resize(quad->weights_.size());
  for (size_t n = 0; n < quad->weights_.size(); ++n)
    quad_data[n] = {quad->weights_[n], quad->abscissae_[n].theta, quad->abscissae_[n].phi};
  return LuaReturn(L, quad_data);
}

} // namespace opensnlua
