// SPDX-FileCopyrightText: 2024 The OpenSn Authors <https://open-sn.github.io/opensn/>
// SPDX-License-Identifier: MIT

#include "lbs_common_lua_functions.h"
#include "modules/linear_boltzmann_solvers/lbs_solver/lbs_solver.h"
#include "framework/runtime.h"
#include "framework/lua.h"
#include "lua/framework/console/console.h"

using namespace opensn;

namespace opensnlua::lbs
{

RegisterLuaFunctionNamespace(LBSComputeFissionRate, lbs, ComputeFissionRate);

int
LBSComputeFissionRate(lua_State* L)
{
  const std::string fname = "lbs.ComputeFissionRate";
  LuaCheckArgs<size_t, std::string>(L, fname);

  // Get pointer to solver
  const auto solver_handle = LuaArg<size_t>(L, 1);
  auto& lbs_solver =
    opensn::GetStackItem<opensn::lbs::LBSSolver>(opensn::object_stack, solver_handle, fname);

  const auto nature = LuaArg<std::string>(L, 2);
  const auto& phi = nature == "OLD" ? lbs_solver.PhiOldLocal() : lbs_solver.PhiNewLocal();

  const double fission_rate = lbs_solver.ComputeFissionRate(phi);
  return LuaReturn(L, fission_rate);
}

} // namespace opensnlua::lbs
