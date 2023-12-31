#include "lbs_do_lua_utils.h"

#include "modules/linear_boltzmann_solvers/b_discrete_ordinates_solver/lbs_discrete_ordinates_solver.h"

#include "framework/console/console.h"

#include "framework/runtime.h"

namespace opensnlua::lbs
{

RegisterLuaFunctionAsIs(chiLBSComputeBalance);

int
chiLBSComputeBalance(lua_State* L)
{
  const std::string fname = "chiLBSComputeBalance";
  const int num_args = lua_gettop(L);

  if (num_args != 1) LuaPostArgAmountError(fname, 1, num_args);

  LuaCheckNilValue(fname, L, 1);

  // Get pointer to solver
  const int solver_handle = lua_tonumber(L, 1);

  auto& lbs_solver = opensn::Chi::GetStackItem<opensn::lbs::DiscreteOrdinatesSolver>(
    opensn::Chi::object_stack, solver_handle, fname);

  lbs_solver.ComputeBalance();

  return 0;
}

} // namespace opensnlua::lbs
