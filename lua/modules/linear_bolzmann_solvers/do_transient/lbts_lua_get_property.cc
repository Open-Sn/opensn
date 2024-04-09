// SPDX-FileCopyrightText: 2024 The OpenSn Authors <https://open-sn.github.io/opensn/>
// SPDX-License-Identifier: MIT

#include "modules/linear_boltzmann_solvers/discrete_ordinates_transient_solver/lbts_transient_solver.h"

#if 0
#include "framework/chi_runtime.h"
#include "framework/logging/chi_log.h"

#define PropertyArgCntErr(prop_name)                                                               \
  throw std::logic_error(fname + ": Insufficient amount of arguments, " +                          \
                         std::to_string(num_args) + ", for property " + prop_name)

namespace lbs::lbts_lua_utils
{

//###################################################################
/**Gets a property of a LBS-TransientSolver or any derived class.

\param SolverIndex int Handle to the solver for which the set is to be created.
\param PropertyName string String name for a specific property.


##_

##ProperyName
"TIMESTEP"\n
Returns the current timestep.\n\n

"TIMESTOP"\n
Returns the stop time.\n\n

"INHIBIT_ADVANCE"\n
Returns the flag.\n\n

"TIME"\n
Returns the simulation time before the latest timestep.\n\n


\author Zachary Hardy*/
int
LBTSGetProperty(lua_State* L)
{
  const std::string fname = "LBTSSetProperty";
  const int num_args = lua_gettop(L);
  if (num_args < 2) LuaPostArgAmountError(fname, 2, num_args);

  // Get the solver
  const int solver_handle = lua_tointeger(L, 1);

  auto& solver =
    chi::GetStackItem<lbs::DiscOrdTransientSolver>(chi::object_stack, solver_handle, fname);

  // Get the property
  LuaCheckStringValue(fname, L, 2);
  const std::string property = lua_tostring(L, 2);

  if (property == "TIMESTEP")
  {
    lua_pushnumber(L, solver.dt_);
    return 1;
  }
  else if (property == "TIMESTOP")
  {
    lua_pushnumber(L, solver.transient_options_.t_final);
    return 1;
  }
  else if (property == "INHIBIT_ADVANCE")
  {
    lua_pushboolean(L, solver.transient_options_.inhibit_advance);
    return 1;
  }
  else if (property == "TIME")
  {
    lua_pushnumber(L, solver.time_);
    return 1;
  }
  else
    throw std::logic_error(fname + ": unsupported property name \"" + property + "\".");
}

} // namespace lbs::lbts_lua_utils
#endif
