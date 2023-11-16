#include "lbsadj_lua_utils.h"

#include "modules/linear_boltzmann_solvers/c_discrete_ordinates_adjoint_solver/lbs_adj_solver.h"

#include "framework/runtime.h"

#include "framework/console/console.h"

namespace opensnlua::lbs
{

RegisterLuaFunctionAsIs(chiAdjointSolverComputeInnerProduct);

int
chiAdjointSolverComputeInnerProduct(lua_State* L)
{
  const std::string fname = __FUNCTION__;
  const int num_args = lua_gettop(L);
  if (num_args != 1) LuaPostArgAmountError(fname, 1, num_args);

  LuaCheckNilValue(fname, L, 1);

  const int solver_handle = lua_tointeger(L, 1);

  auto& solver = opensn::Chi::GetStackItem<opensn::lbs::DiscreteOrdinatesAdjointSolver>(
    opensn::Chi::object_stack, solver_handle, fname);

  const double ip_Q_phi_star = solver.ComputeInnerProduct();

  lua_pushnumber(L, ip_Q_phi_star);
  return 1;
}

} // namespace opensnlua::lbs
