#include "modules/LinearBoltzmannSolvers/A_LBSSolver/lbs_solver.h"

#include "framework/physics/FieldFunction/fieldfunction_gridbased.h"

#include "framework/chi_runtime.h"
#include "framework/logging/chi_log.h"

namespace lbs::common_lua_utils
{

/**Obtains a list of field functions, related only to scalar flux,
from the transport solver.

\param SolverIndex int Handle to the solver for which the list is to be
obtained.

\return Pair Table and count. Returns an array of handles and the amount of
elements in it (indexed from 1). \ingroup LBSLuaFunctions \author Jan*/
int
chiLBSGetScalarFieldFunctionList(lua_State* L)
{
  const std::string fname = __FUNCTION__;
  const int num_args = lua_gettop(L);
  if (num_args != 1) LuaPostArgAmountError(fname, 1, num_args);

  LuaCheckNilValue(fname, L, 1);
  LuaCheckNumberValue(fname, L, 1);

  //============================================= Get pointer to solver
  const int solver_handle = lua_tonumber(L, 1);
  const auto& lbs_solver =
    Chi::GetStackItem<lbs::LBSSolver>(Chi::object_stack, solver_handle, fname);

  /**Lambda for matching a field function smart pointer to one on
   * the runtime stack.*/
  auto GetStackFFHandle = [](std::shared_ptr<chi_physics::FieldFunctionGridBased>& local_ff)
  {
    size_t stack_ff_counter = 0;
    for (auto& stack_ff : Chi::field_function_stack)
    {
      if (stack_ff == local_ff) return stack_ff_counter;

      ++stack_ff_counter;
    }

    ChiLogicalError("Scalar field function lookup error");
  };

  //======================================== Building table of handles
  lua_newtable(L);
  lua_Integer count = 0;

  // Flux moments first
  for (int g = 0; g < lbs_solver.NumGroups(); g++)
  {
    for (int m = 0; m < lbs_solver.NumMoments(); m++)
    {
      const size_t ff = lbs_solver.MapPhiFieldFunction(g, m);
      auto local_ff = lbs_solver.GetFieldFunctions()[ff];

      if (m != 0) continue;

      lua_pushinteger(L, 1 + count++);
      lua_pushinteger(L, static_cast<lua_Integer>(GetStackFFHandle(local_ff)));

      lua_settable(L, -3);
    }
  }

  //// Power generation
  // if (lbs_solver.Options().power_field_function_on)
  //{
  //   const size_t ff = lbs_solver.GetHandleToPowerGenFieldFunc();
  //   auto local_ff = lbs_solver.GetFieldFunctions()[ff];
  //
  //   lua_pushinteger(L, 1 + count++);
  //   lua_pushinteger(L, static_cast<lua_Integer>(GetStackFFHandle(local_ff)));
  //
  //   lua_settable(L, -3);
  // }

  lua_pushinteger(L, count);
  return 2;
}

} // namespace lbs::common_lua_utils
