// SPDX-FileCopyrightText: 2024 The OpenSn Authors <https://open-sn.github.io/opensn/>
// SPDX-License-Identifier: MIT

#include "lua/framework/lua.h"
#include "lua/framework/console/console.h"
#include "lua/framework/math/quadratures/sldfesq/sldfe.h"
#include "framework/math/quadratures/angular/sldfe_sq_quadrature.h"
#include "framework/logging/log.h"
#include "framework/runtime.h"

using namespace opensn;

namespace opensnlua
{

RegisterLuaFunctionInNamespace(LocallyRefineSLDFESQAngularQuadrature, aquad, LocallyRefineSLDFESQ);

int
LocallyRefineSLDFESQAngularQuadrature(lua_State* L)
{
  LuaCheckArgs<size_t, Vector3, double>(L, "aquad.LocallyRefineSLDFESQ");

  auto handle = LuaArg<size_t>(L, 1);
  auto ref_dir = LuaArg<Vector3>(L, 2);
  auto cone_size = LuaArg<double>(L, 3);
  bool ref_dir_as_plane_normal = LuaArgOptional<bool>(L, 4, false);

  try
  {
    auto ref_quadrature = opensn::angular_quadrature_stack.at(handle);
    if (ref_quadrature->type_ == AngularQuadratureType::SLDFESQ)
    {
      auto sldfesq = std::dynamic_pointer_cast<SimplifiedLDFESQ::Quadrature>(ref_quadrature);

      sldfesq->LocallyRefine(ref_dir, cone_size, ref_dir_as_plane_normal);
    }
    else
    {
      opensn::log.LogAllError() << "LocallyRefineSLDFESQAngularQuadrature: "
                                   "Invalid angular quadrature type.";
      opensn::Exit(EXIT_FAILURE);
    }
  }
  catch (const std::out_of_range& o)
  {
    opensn::log.LogAllError() << "LocallyRefineSLDFESQAngularQuadrature: "
                                 "Invalid handle to angular quadrature.";
    opensn::Exit(EXIT_FAILURE);
  }
  catch (...)
  {
    opensn::log.LogAllError() << "LocallyRefineSLDFESQAngularQuadrature: "
                                 "Call failed with unknown error.";
    opensn::Exit(EXIT_FAILURE);
  }

  return LuaReturn(L);
}

} // namespace opensnlua
