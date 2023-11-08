#include "lbs_DO_lua_utils.h"

#include "modules/LinearBoltzmannSolvers/B_DiscreteOrdinatesSolver/lbs_discrete_ordinates_solver.h"

#include "framework/console/chi_console.h"

#include "framework/chi_runtime.h"

namespace lbs::disc_ord_lua_utils
{

RegisterLuaFunctionAsIs(chiLBSComputeBalance);

// ###################################################################
/**Computes balance tables and prints it to the console.
 *
\param SolverIndex int Handle to the solver for which the list is to be
obtained.

\ingroup LBSLuaFunctions
\author Jan*/
int
chiLBSComputeBalance(lua_State* L)
{
  const std::string fname = "chiLBSComputeBalance";
  const int num_args = lua_gettop(L);

  if (num_args != 1) LuaPostArgAmountError(fname, 1, num_args);

  LuaCheckNilValue(fname, L, 1);

  //============================================= Get pointer to solver
  const int solver_handle = lua_tonumber(L, 1);

  auto& lbs_solver =
    Chi::GetStackItem<lbs::DiscreteOrdinatesSolver>(Chi::object_stack, solver_handle, fname);

  lbs_solver.ComputeBalance();

  return 0;
}

} // namespace lbs::disc_ord_lua_utils
