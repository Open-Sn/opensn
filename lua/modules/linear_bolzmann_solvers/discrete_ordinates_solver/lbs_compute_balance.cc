#include "lbs_do_lua_utils.h"

#include "modules/linear_boltzmann_solvers/discrete_ordinates_solver/lbs_discrete_ordinates_solver.h"

#include "framework/console/console.h"

#include "framework/runtime.h"

namespace opensnlua::lbs
{

RegisterLuaFunctionNamespace(LBSComputeBalance, lbs, ComputeBalance);

int
LBSComputeBalance(lua_State* L)
{
  const std::string fname = "LBSComputeBalance";
  const int num_args = lua_gettop(L);

  if (num_args != 1)
    LuaPostArgAmountError(fname, 1, num_args);

  // Get pointer to solver
  const auto solver_handle = LuaArg<size_t>(L, 1);

  auto& lbs_solver = opensn::GetStackItem<opensn::lbs::DiscreteOrdinatesSolver>(
    opensn::object_stack, solver_handle, fname);

  lbs_solver.ComputeBalance();

  return LuaReturn(L);
}

} // namespace opensnlua::lbs
