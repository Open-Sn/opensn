// SPDX-FileCopyrightText: 2024 The OpenSn Authors <https://open-sn.github.io/opensn/>
// SPDX-License-Identifier: MIT

#include "modules/linear_boltzmann_solvers/discrete_ordinates_transient_solver/lbts_transient_solver.h"
#if 0
#include "framework/runtime.h"
#include "framework/logging/log.h"

#define PropertyArgCntErr(prop_name)                                                               \
  throw std::logic_error(fname + ": Insufficient amount of arguments, " +                          \
                         std::to_string(num_args) + ", for property " + prop_name)

namespace lbts_lua_utils
{

//###################################################################
/**Sets a property of the LBS-TransientSolver or any derived class.

\param SolverIndex int Handle to the solver for which the set is to be created.
\param PropertyName string String name for a specific property.


##_

##ProperyName
"TIMESTEP"\n
Expects to be followed by a floating point value representing the timestep.\n\n

"TIME"\n
Sets the current time. [Default=0.0]\n\n

"TIMESTOP"\n
Sets the termination time. [Default=0.1]\n\n

"MAX_TIMESTEPS"\n
Sets the maximum number of timesteps during a call to Execute. A negative
value disable the logic to check for timestep-counters exceeding this amount.
[Default=10]\n\n

"INHIBIT_ADVANCE"\n
Sets a flag that prevents the time from advancing. Can be used to reject a
timestep. [Default=false]\n\n

"VERBOSITY_LEVEL"\n
Sets the verbosity level. Level 0, nothing gets printed. Level 1, basic info.
 [Default=1]\n\n

"TIMESTEP_METHOD"\n
Sets the time-stepping method. Can be "CRANK_NICHOLSON" or "BACKWARD_EULER".
[Default="CRANK_NICHOLSON"]\n\n

"CALLBACK"\n
Sets the timestep callback function. [Default=Nothing]\n\n

"SCALE_FISSION_XS"\n
Sets a boolean flag to normalize the fission cross sections by the
k-eigenvalue. [Default=false]\n\n

"NORMALIZATION_METHOD"\n
Sets the initial data normalization data. Can be "TOTAL_POWER",
"POWER_DENSITY", or "NONE". [Default="TOTAL_POWER"]\n\n

\author Zachary Hardy*/
int
LBTSSetProperty(lua_State* L)
{
  const std::string fname = "LBTSSetProperty";
  const int num_args = lua_gettop(L);
  if (num_args < 2) LuaPostArgAmountError(fname, 2, num_args);

  // Get the solver
  LuaCheckNilValue(fname, L, 1);
  const int solver_handle = lua_tointeger(L, 1);

  auto& solver =
    opensn::GetStackItem<DiscOrdTransientSolver>(opensn::object_stack, solver_handle, fname);

  // Get the property
  LuaCheckStringValue(fname, L, 2);
  const std::string property = lua_tostring(L, 2);

  if (property == "TIMESTEP")
  {
    if (num_args != 3) PropertyArgCntErr("TIMESTEP");

    LuaCheckNilValue(fname, L, 3);

    const double dt_input = lua_tonumber(L, 3);

    solver.dt_ = dt_input;

    opensn::log.Log() << solver.TextName() << ": dt set to " << std::to_string(dt_input);
  }
  else if (property == "TIME")
  {
    if (num_args != 3) PropertyArgCntErr("TIME");

    LuaCheckNilValue(fname, L, 3);

    const double t_input = lua_tonumber(L, 3);

    solver.time_ = t_input;

    opensn::log.Log() << solver.TextName() << ": time set to " << std::to_string(t_input);
  }
  else if (property == "TIMESTOP")
  {
    if (num_args != 3) PropertyArgCntErr("TIMESTOP");

    LuaCheckNilValue(fname, L, 3);

    const double t_input = lua_tonumber(L, 3);

    solver.transient_options_.t_final = t_input;

    opensn::log.Log() << solver.TextName() << ": t_final set to " << std::to_string(t_input);
  }
  else if (property == "MAX_TIMESTEPS")
  {
    if (num_args != 3) PropertyArgCntErr("MAX_TIMESTEPS");

    LuaCheckNilValue(fname, L, 3);

    const int t_input = lua_tointeger(L, 3);

    solver.transient_options_.max_time_steps = t_input;

    opensn::log.Log() << solver.TextName() << ": max_time_steps set to " << std::to_string(t_input);
  }
  else if (property == "INHIBIT_ADVANCE")
  {
    if (num_args != 3) PropertyArgCntErr("INHIBIT_ADVANCE");

    LuaCheckNilValue(fname, L, 3);

    const bool inhibit_advance = lua_toboolean(L, 3);

    solver.transient_options_.inhibit_advance = inhibit_advance;

    opensn::log.Log() << solver.TextName() << ": inhibit_advance set to "
                   << std::to_string(inhibit_advance);
  }
  else if (property == "VERBOSITY_LEVEL")
  {
    if (num_args != 3) PropertyArgCntErr("VERBOSITY_LEVEL");

    LuaCheckNilValue(fname, L, 3);

    const int verbosity_level = lua_tointeger(L, 3);

    solver.transient_options_.verbosity_level = verbosity_level;

    opensn::log.Log() << solver.TextName() << ": verbosity_level set to "
                   << std::to_string(verbosity_level);
  }
  else if (property == "TIMESTEP_METHOD")
  {
    if (num_args != 3) PropertyArgCntErr("TIMESTEP_METHOD");

    LuaCheckNilValue(fname, L, 3);

    const std::string option = lua_tostring(L, 3);

    if (option == "BACKWARD_EULER") solver.method = opensn::SteppingMethod::IMPLICIT_EULER;
    else if (option == "CRANK_NICHOLSON")
      solver.method = opensn::SteppingMethod::CRANK_NICOLSON;
    else
      throw std::invalid_argument(fname +
                                  ": Only the following timestepping "
                                  "methods are supported: \"CRANK_NICHOLSON\", \"BACKWARD_EULER\"");

    opensn::log.Log() << solver.TextName() << ": method set to " << option;
  }
  else if (property == "CALLBACK")
  {
    if (num_args != 3) PropertyArgCntErr("CALLBACK");

    LuaCheckNilValue(fname, L, 3);

    const std::string cbfname = lua_tostring(L, 3);

    solver.transient_options_.console_call_back_function = cbfname;
    opensn::log.Log() << solver.TextName() << ": console_call_back_function set to " << cbfname;
  }
  else if (property == "SCALE_FISSION_XS")
  {
    if (num_args != 3) PropertyArgCntErr("SCALE_FISSION_XS");

    LuaCheckNilValue(fname, L, 3);

    const bool scale_fission_xs = lua_toboolean(L, 3);
    solver.transient_options_.scale_fission_xs = scale_fission_xs;

    opensn::log.Log() << solver.TextName() << ": scale_fission_xs set to "
                   << std::to_string(scale_fission_xs);
  }
  else if (property == "NORMALIZATION_METHOD")
  {
    if (num_args != 3) PropertyArgCntErr("NORMALIZATION_METHOD");

    LuaCheckNilValue(fname, L, 3);

    const std::string option = lua_tostring(L, 3);

    if (option == "TOTAL_POWER")
      solver.transient_options_.normalization_method =
        DiscOrdTransientSolver::NormalizationMethod::TOTAL_POWER;
    else if (option == "POWER_DENSITY")
      solver.transient_options_.normalization_method =
        DiscOrdTransientSolver::NormalizationMethod::POWER_DENSITY;
    else if (option == "NONE")
      solver.transient_options_.normalization_method =
        DiscOrdTransientSolver::NormalizationMethod::NONE;
    else
      throw std::invalid_argument(fname + ": Only the following normalization methods are " +
                                  "supported: \"TOTAL_POWER\", \"POWER_DENSITY\", \"NONE\"");

    opensn::log.Log() << solver.TextName() << ": normalization_method set to " << option;
  }
  else
    throw std::logic_error(fname + ": unsupported property name \"" + property + "\".");

  return 0;
}

} // namespace lbts_lua_utils
#endif
