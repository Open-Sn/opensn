// SPDX-FileCopyrightText: 2024 The OpenSn Authors <https://open-sn.github.io/opensn/>
// SPDX-License-Identifier: MIT

#include "lua/modules/linear_bolzmann_solvers/lbs_solver/lbs_common_lua_functions.h"
#include "lua/framework/lua.h"
#include "lua/framework/console/console.h"
#include "modules/linear_boltzmann_solvers/lbs_solver/lbs_solver.h"
#include "framework/field_functions/field_function_grid_based.h"
#include "framework/runtime.h"
#include "framework/logging/log.h"

using namespace opensn;

namespace opensnlua
{

RegisterLuaFunctionInNamespace(LBSGetScalarFieldFunctionList, lbs, GetScalarFieldFunctionList);
RegisterLuaFunctionInNamespace(LBSGetPowerFieldFunction, lbs, GetPowerFieldFunction);

int
LBSGetScalarFieldFunctionList(lua_State* L)
{
  const std::string fname = "lbs.GetScalarFieldFunctionList";
  LuaCheckArgs<size_t>(L, fname);

  // Get pointer to solver
  const auto solver_handle = LuaArg<size_t>(L, 1);
  const auto& lbs_solver =
    opensn::GetStackItem<opensn::LBSSolver>(opensn::object_stack, solver_handle, fname);

  /**Lambda for matching a field function smart pointer to one on
   * the runtime stack.*/
  auto GetStackFFHandle = [](std::shared_ptr<FieldFunctionGridBased>& local_ff)
  {
    size_t stack_ff_counter = 0;
    for (auto& stack_ff : opensn::field_function_stack)
    {
      if (stack_ff == local_ff)
        return stack_ff_counter;

      ++stack_ff_counter;
    }

    OpenSnLogicalError("Scalar field function lookup error");
  };

  // Building table of handles
  std::vector<size_t> ff_handles;
  // Flux moments first
  for (int g = 0; g < lbs_solver.NumGroups(); g++)
  {
    for (int m = 0; m < lbs_solver.NumMoments(); m++)
    {
      const size_t ff = lbs_solver.MapPhiFieldFunction(g, m);
      auto local_ff = lbs_solver.GetFieldFunctions()[ff];

      if (m != 0)
        continue;

      ff_handles.push_back(GetStackFFHandle(local_ff));
    }
  }

  return LuaReturn(L, ff_handles, ff_handles.size());
}

int
LBSGetPowerFieldFunction(lua_State* L)
{
  const std::string fname = "lbs.GetPowerFieldFunction";
  LuaCheckArgs<size_t>(L, fname);

  // Get pointer to solver
  const auto solver_handle = LuaArg<size_t>(L, 1);
  const auto& lbs_solver =
    opensn::GetStackItem<opensn::LBSSolver>(opensn::object_stack, solver_handle, fname);

  if (lbs_solver.Options().power_field_function_on)
    return LuaReturn(L, lbs_solver.GetHandleToPowerGenFieldFunc());
  else
    throw std::logic_error("The power field function is not enabled. Use \"power_field_function_on "
                           "= true\" in the input options to enable it.");
}
} // namespace opensnlua
