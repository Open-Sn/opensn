#include "opensn/modules/LinearBoltzmannSolvers/A_LBSSolver/lbs_solver.h"

#include "opensn/framework/chi_runtime.h"

namespace lbs::common_lua_utils
{

//###################################################################
/**Computes and returns the fission rate.
 *
\param SolverIndex int Handle to the solver maintaining the information.
\param OldNewOption string "NEW" or "OLD". For steady state solvers, the
                           "OLD" option would give the fission rate for
                           the previous iterate. [Default="NEW"]

\return double The fission rate.

\ingroup LBSLuaFunctions
\author Jan*/
int
chiLBSComputeFissionRate(lua_State* L)
{
  const std::string fname = "chiLBSComputeFissionRate";
  const int num_args = lua_gettop(L);

  if (num_args != 2) LuaPostArgAmountError(fname, 2, num_args);

  LuaCheckNilValue(fname, L, 1);

  //============================================= Get pointer to solver
  const int solver_handle = lua_tonumber(L, 1);

  auto& lbs_solver = Chi::GetStackItem<lbs::LBSSolver>(Chi::object_stack, solver_handle, fname);

  LuaCheckStringValue(fname, L, 2);
  const std::string nature = lua_tostring(L, 2);
  const auto& phi = nature == "OLD" ? lbs_solver.PhiOldLocal() : lbs_solver.PhiNewLocal();

  const double fission_rate = lbs_solver.ComputeFissionRate(phi);

  lua_pushnumber(L, fission_rate);

  return 1;
}

} // namespace lbs::common_lua_utils
