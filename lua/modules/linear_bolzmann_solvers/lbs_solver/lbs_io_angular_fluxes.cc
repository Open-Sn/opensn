// SPDX-FileCopyrightText: 2024 The OpenSn Authors <https://open-sn.github.io/opensn/>
// SPDX-License-Identifier: MIT

#include "lua/modules/linear_bolzmann_solvers/lbs_solver/lbs_common_lua_functions.h"
#include "lua/framework/lua.h"
#include "lua/framework/console/console.h"
#include "modules/linear_boltzmann_solvers/lbs_solver/lbs_solver.h"
#include "modules/linear_boltzmann_solvers/lbs_solver/groupset/lbs_groupset.h"
#include "modules/linear_boltzmann_solvers/lbs_solver/io/lbs_solver_io.h"
#include "framework/runtime.h"
#include "framework/logging/log.h"

using namespace opensn;

namespace opensnlua
{

RegisterLuaFunctionInNamespace(LBSWriteGroupsetAngularFlux, lbs, WriteGroupsetAngularFlux);
RegisterLuaFunctionInNamespace(LBSReadGroupsetAngularFlux, lbs, ReadGroupsetAngularFlux);

RegisterLuaFunctionInNamespace(LBSWriteAngularFluxes, lbs, WriteAngularFluxes);
RegisterLuaFunctionInNamespace(LBSReadAngularFluxes, lbs, ReadAngularFluxes);

int
LBSWriteGroupsetAngularFlux(lua_State* L)
{
  const std::string fname = "lbs.WriteGroupsetAngularFlux";
  LuaCheckArgs<size_t, int, std::string>(L, fname);

  const auto solver_handle = LuaArg<size_t>(L, 1);
  const auto groupset_index = LuaArg<int>(L, 2);
  const auto file_base = LuaArg<std::string>(L, 3);

  // Get pointer to solver
  auto& lbs_solver =
    opensn::GetStackItem<opensn::LBSSolver>(opensn::object_stack, solver_handle, fname);
  LBSSolverIO::WriteGroupsetAngularFluxes(lbs_solver, groupset_index, file_base);

  return LuaReturn(L);
}

int
LBSReadGroupsetAngularFlux(lua_State* L)
{
  const std::string fname = "lbs.ReadGroupsetAngularFlux";
  LuaCheckArgs<size_t, int, std::string>(L, fname);

  const auto solver_handle = LuaArg<size_t>(L, 1);
  const auto groupset_index = LuaArg<int>(L, 2);
  const auto file_base = LuaArg<std::string>(L, 3);

  // Get pointer to solver
  auto& lbs_solver =
    opensn::GetStackItem<opensn::LBSSolver>(opensn::object_stack, solver_handle, fname);
  LBSSolverIO::ReadGroupsetAngularFluxes(lbs_solver, groupset_index, file_base);

  return LuaReturn(L);
}

int
LBSWriteAngularFluxes(lua_State* L)
{
  const std::string fname = "lbs.WriteAngularFluxes";
  LuaCheckArgs<size_t, std::string>(L, fname);

  const auto solver_handle = LuaArg<size_t>(L, 1);
  const auto file_base = LuaArg<std::string>(L, 2);

  // Get pointer to solver
  auto& lbs_solver =
    opensn::GetStackItem<opensn::LBSSolver>(opensn::object_stack, solver_handle, fname);
  LBSSolverIO::WriteAngularFluxes(lbs_solver, file_base);

  return LuaReturn(L);
}

int
LBSReadAngularFluxes(lua_State* L)
{
  const std::string fname = "lbs.ReadAngularFluxes";
  LuaCheckArgs<size_t, std::string>(L, fname);

  const auto solver_handle = LuaArg<size_t>(L, 1);
  const auto file_base = LuaArg<std::string>(L, 2);

  // Get pointer to solver
  auto& lbs_solver =
    opensn::GetStackItem<opensn::LBSSolver>(opensn::object_stack, solver_handle, fname);
  LBSSolverIO::ReadAngularFluxes(lbs_solver, file_base);

  return LuaReturn(L);
}

} // namespace opensnlua
