// SPDX-FileCopyrightText: 2024 The OpenSn Authors <https://open-sn.github.io/opensn/>
// SPDX-License-Identifier: MIT

#include "lua/modules/linear_bolzmann_solvers/discrete_ordinates_adjoint_solver/lbsadj_lua_utils.h"
#include "modules/linear_boltzmann_solvers/discrete_ordinates_adjoint_solver/lbs_adj_solver.h"
#include "framework/runtime.h"
#include "lua/framework/console/console.h"

using namespace opensn;

namespace opensnlua::lbs
{

RegisterLuaFunction(AdjointSolverExportImportanceMapBinary);

int
AdjointSolverExportImportanceMapBinary(lua_State* L)
{
  const std::string fname = "AdjointSolverExportImportanceMapBinary";
  LuaCheckArgs<int, std::string>(L, fname);

  const auto solver_handle = LuaArg<int>(L, 1);
  const auto file_name = LuaArg<std::string>(L, 2);

  auto& solver = opensn::GetStackItem<opensn::lbs::DiscreteOrdinatesAdjointSolver>(
    opensn::object_stack, solver_handle, fname);

  solver.ExportImportanceMap(file_name);

  return LuaReturn(L);
}

} // namespace opensnlua::lbs
