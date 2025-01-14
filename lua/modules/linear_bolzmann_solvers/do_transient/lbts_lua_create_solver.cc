// SPDX-FileCopyrightText: 2024 The OpenSn Authors <https://open-sn.github.io/opensn/>
// SPDX-License-Identifier: MIT

#include "modules/linear_boltzmann_solvers/discrete_ordinates_transient_solver/lbts_transient_solver.h"
#include "framework/runtime.h"

#if 0
namespace lbts_lua_utils
{

//###################################################################
/**Creates a LBS-TransientSolver.

\param SolverName string Optional. The text name of the solver.
                         [Default="LBSTransientSolver"].

*/
int
LBSCreateTransientSolver(lua_State* L)
{
  const std::string fname = "LBSCreateTransientSolver";
  const int num_args = lua_gettop(L);

  std::string solver_name = "LBSTransient";
  if (num_args == 1)
  {
    LuaCheckStringValue(fname, L, 1);
    solver_name = lua_tostring(L, 1);
  }

  auto new_solver = std::make_shared<DiscOrdTransientSolver>(solver_name);

  opensn::object_stack.push_back(new_solver);

  lua_pushinteger(L, static_cast<lua_Integer>(opensn::object_stack.size() - 1));
  return 1;
}

} // namespace lbts_lua_utils
#endif
