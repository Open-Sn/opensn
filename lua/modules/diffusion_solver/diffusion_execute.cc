#include "opensn/framework/chi_lua.h"
#include "opensn/modules/DiffusionSolver/Solver/diffusion_solver.h"

#include "opensn/framework/chi_runtime.h"

//#############################################################################
/** Initialize the Diffusion solver.
 *
\param SolverHandle int Handle to an existing diffusion solver.

\return Success bool Returns if initialization failed.
\ingroup LuaDiffusion
\author Jan*/
int
chiDiffusionExecute(lua_State* L)
{
  const size_t solver_index = lua_tonumber(L, 1);
  auto& solver =
    Chi::GetStackItem<chi_diffusion::Solver>(Chi::object_stack, solver_index, __FUNCTION__);

  solver.ExecuteS();

  return 0;
}
