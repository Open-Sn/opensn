#include "modules/LinearBoltzmannSolvers/D_DO_Transient/lbts_transient_solver.h"
#include "lbts_lua_utils.h"
#include "framework/runtime.h"

#if 0
namespace lbs::lbts_lua_utils
{

//###################################################################
/**Creates a LBS-TransientSolver.

\param SolverName string Optional. The text name of the solver.
                         [Default="LBSTransientSolver"].

\author Zachary Hardy*/
int
chiLBSCreateTransientSolver(lua_State* L)
{
  const std::string fname = "chiLBSCreateTransientSolver";
  const int num_args = lua_gettop(L);

  std::string solver_name = "LBSTransient";
  if (num_args == 1)
  {
    LuaCheckStringValue(fname, L, 1);
    solver_name = lua_tostring(L, 1);
  }

  auto new_solver = std::make_shared<lbs::DiscOrdTransientSolver>(solver_name);

  chi::object_stack.push_back(new_solver);

  lua_pushinteger(L, static_cast<lua_Integer>(chi::object_stack.size() - 1));
  return 1;
}

} // namespace lbs::lbts_lua_utils
#endif
