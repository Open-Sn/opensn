// SPDX-FileCopyrightText: 2024 The OpenSn Authors <https://open-sn.github.io/opensn/>
// SPDX-License-Identifier: MIT

#include "framework/lua.h"
#include "framework/math/quadratures/angular/sldfe_sq_quadrature.h"
#include "framework/console/console.h"
#include "framework/runtime.h"
#include "sldfe_lua.h"

using namespace opensn;

namespace opensnlua
{

RegisterLuaFunctionInNamespace(CreateSLDFESQAngularQuadrature,
                               aquad,
                               CreateSLDFESQAngularQuadrature);

int
CreateSLDFESQAngularQuadrature(lua_State* L)
{
  LuaCheckArgs<int>(L, "aquad.CreateSLDFESQAngularQuadrature");

  auto init_refinement_level = LuaArg<int>(L, 1);

  auto sldfesq = new SimplifiedLDFESQ::Quadrature;
  sldfesq->GenerateInitialRefinement(init_refinement_level);

  std::shared_ptr<AngularQuadrature> new_ang_quad =
    std::shared_ptr<SimplifiedLDFESQ::Quadrature>(sldfesq);

  opensn::angular_quadrature_stack.push_back(new_ang_quad);
  const size_t index = opensn::angular_quadrature_stack.size() - 1;
  return LuaReturn(L, index);
}

} // namespace opensnlua
