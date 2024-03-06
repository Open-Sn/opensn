#include "lbsadj_lua_utils.h"

#include "modules/linear_boltzmann_solvers/discrete_ordinates_adjoint_solver/lbs_adj_solver.h"

#include "framework/runtime.h"

#include "framework/console/console.h"

namespace opensnlua::lbs
{

RegisterLuaFunctionAsIs(AdjointSolverComputeInnerProduct);

int
AdjointSolverComputeInnerProduct(lua_State* L)
{
  const std::string fname = __FUNCTION__;
  const int num_args = lua_gettop(L);
  if (num_args != 1)
    LuaPostArgAmountError(fname, 1, num_args);

  LuaCheckNilValue(fname, L, 1);

  const int solver_handle = lua_tointeger(L, 1);
  auto& solver = opensn::GetStackItem<opensn::lbs::DiscreteOrdinatesAdjointSolver>(
    opensn::object_stack, solver_handle, fname);

  const double response = solver.ComputeInnerProduct();
  lua_pushnumber(L, response);
  return 1;
}

} // namespace opensnlua::lbs
