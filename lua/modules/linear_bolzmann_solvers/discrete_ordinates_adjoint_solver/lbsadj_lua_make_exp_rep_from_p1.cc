// SPDX-FileCopyrightText: 2024 The OpenSn Authors <https://open-sn.github.io/opensn/>
// SPDX-License-Identifier: MIT

#include "lua/modules/linear_bolzmann_solvers/discrete_ordinates_adjoint_solver/lbsadj_lua_utils.h"
#include "modules/linear_boltzmann_solvers/discrete_ordinates_adjoint_solver/lbs_adjoint.h"
#include "framework/console/console.h"
#include <stdexcept>

using namespace opensn;

namespace opensnlua::lbs
{

RegisterLuaFunction(AdjointSolverMakeExpRepFromP1Moments);

int
AdjointSolverMakeExpRepFromP1Moments(lua_State* L)
{
  const std::string fname = "AdjointSolverMakeExpRepFromP1Moments";
  LuaCheckArgs<std::vector<double>>(L, fname);

  auto P1 = LuaArg<std::vector<double>>(L, 1);
  if (P1.size() != 4)
    throw std::invalid_argument(fname + ": Supplied table argument must have 4 entries.");

  auto verbose = LuaArgOptional<bool>(L, 2, false);

  auto solution = opensn::lbs::MakeExpRepFromP1({P1[0], P1[1], P1[2], P1[3]}, verbose);

  return LuaReturn(L, solution[0], solution[1]);
}
} // namespace opensnlua::lbs
