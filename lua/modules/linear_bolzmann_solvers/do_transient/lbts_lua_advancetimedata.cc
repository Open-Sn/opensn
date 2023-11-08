#include "modules/LinearBoltzmannSolvers/D_DO_Transient/lbts_transient_solver.h"

#include "framework/chi_runtime.h"
#include "framework/logging/chi_log.h"

#if 0
#define PropertyArgCntErr(prop_name)                                                               \
  throw std::logic_error(fname + ": Insufficient amount of arguments, " +                          \
                         std::to_string(num_args) + ", for property " + prop_name)

namespace lbs::lbts_lua_utils
{

//###################################################################
/**Advances time dependent data to the next timestep.

\param SolverIndex int Handle to the solver..


\author Zachary Hardy*/
int
chiLBTSAdvanceTimeData(lua_State* L)
{
  const std::string fname = "chiLBTSSetProperty";
  const int num_args = lua_gettop(L);
  if (num_args != 1) LuaPostArgAmountError(fname, 1, num_args);

  // Get the solver
  LuaCheckNilValue(fname, L, 1);
  const int solver_handle = lua_tointeger(L, 1);

  auto& solver =
    chi::GetStackItem<lbs::DiscOrdTransientSolver>(chi::object_stack, solver_handle, fname);

  solver.Advance();

  return 0;
}

} // namespace lbs::lbts_lua_utils
#endif
