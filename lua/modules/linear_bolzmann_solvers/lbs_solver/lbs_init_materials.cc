// SPDX-FileCopyrightText: 2024 The OpenSn Authors <https://open-sn.github.io/opensn/>
// SPDX-License-Identifier: MIT

#include "lua/modules/linear_bolzmann_solvers/lbs_solver/lbs_common_lua_functions.h"
#include "modules/linear_boltzmann_solvers/lbs_solver/lbs_solver.h"
#include "framework/runtime.h"
#include "lua/framework/lua.h"
#include "lua/framework/console/console.h"

using namespace opensn;

namespace opensnlua::lbs
{

RegisterLuaFunctionInNamespace(LBSInitializeMaterials, lbs, InitializeMaterials);

int
LBSInitializeMaterials(lua_State* L)
{
  const std::string fname = "lbs.InitializeMaterials";
  LuaCheckArgs<size_t>(L, fname);

  // Get pointer to solver
  const auto solver_handle = LuaArg<size_t>(L, 1);
  auto& lbs_solver =
    opensn::GetStackItem<opensn::lbs::LBSSolver>(opensn::object_stack, solver_handle, fname);

  lbs_solver.InitializeMaterials();

  return LuaReturn(L);
}

} // namespace opensnlua::lbs
