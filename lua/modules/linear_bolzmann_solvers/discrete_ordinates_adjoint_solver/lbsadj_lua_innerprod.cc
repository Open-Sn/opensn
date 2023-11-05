#include "lbsadj_lua_utils.h"

#include "C_DiscreteOrdinatesAdjointSolver/lbsadj_solver.h"

#include "chi_runtime.h"

#include "console/chi_console.h"

namespace lbs::adjoint_lua_utils
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

  auto& solver =
    Chi::GetStackItem<lbs::DiscreteOrdinatesAdjointSolver>(Chi::object_stack, solver_handle, fname);

  const double ip_Q_phi_star = solver.ComputeInnerProduct();

  lua_pushnumber(L, ip_Q_phi_star);
  return 1;
}

} // namespace lbs::adjoint_lua_utils