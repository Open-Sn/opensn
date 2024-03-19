#include "lbsadj_lua_utils.h"

#include "modules/linear_boltzmann_solvers/discrete_ordinates_adjoint_solver/lbs_adjoint.h"

#include <stdexcept>

#include "framework/console/console.h"

using namespace opensn;

namespace opensnlua::lbs
{

RegisterLuaFunctionAsIs(AdjointSolverMakeExpRepFromP1Moments);

int
AdjointSolverMakeExpRepFromP1Moments(lua_State* L)
{
  const std::string fname = __FUNCTION__;
  const int num_args = lua_gettop(L);

  if (num_args < 1)
    LuaPostArgAmountError(fname, 1, num_args);

  LuaCheckNilValue(fname, L, 1);
  LuaCheckTableValue(fname, L, 1);

  std::vector<double> P1;
  LuaPopulateVectorFrom1DArray(fname, L, 1, P1);

  if (P1.size() != 4)
    throw std::invalid_argument(fname + ": Supplied table argument must"
                                        " have 4 entries.");

  auto verbose = LuaArgOptional<bool>(L, 2, false);

  auto solution = opensn::lbs::MakeExpRepFromP1({P1[0], P1[1], P1[2], P1[3]}, verbose);

  lua_pushnumber(L, solution[0]);
  lua_pushnumber(L, solution[1]);
  return 2;
}
} // namespace opensnlua::lbs
