// SPDX-FileCopyrightText: 2024 The OpenSn Authors <https://open-sn.github.io/opensn/>
// SPDX-License-Identifier: MIT

#include "lua/framework/lua.h"
#include "lua/framework/console/console.h"
#include "lua/framework/math/quadratures/quadratures.h"
#include "framework/math/quadratures/angular/angular_quadrature.h"
#include "framework/logging/log.h"
#include "framework/runtime.h"

using namespace opensn;

namespace opensnlua
{

RegisterLuaFunctionInNamespace(OptimizeAngularQuadratureForPolarSymmetry,
                               aquad,
                               OptimizeForPolarSymmetry);

int
OptimizeAngularQuadratureForPolarSymmetry(lua_State* L)
{
  const std::string fname = "aquad.OptimizeForPolarSymmetry";
  LuaCheckArgs<int>(L, fname);

  const auto handle = LuaArg<int>(L, 1);
  auto normalization = LuaArgOptional<double>(L, 2, -1.0);

  auto& quadrature =
    opensn::GetStackItem<AngularQuadrature>(opensn::angular_quadrature_stack, handle, fname);

  if (normalization > 0.0)
    opensn::log.Log() << "Optimizing angular quadrature for polar symmetry. using "
                      << "normalization factor " << normalization << ".";

  quadrature.OptimizeForPolarSymmetry(normalization);

  return LuaReturn(L);
}

} // namespace opensnlua
