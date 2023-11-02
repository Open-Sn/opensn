#include "chi_lua.h"

#include "../Solver/diffusion_solver.h"

#include "chi_runtime.h"

//#############################################################################
/** Initialize the Diffusion solver.
 *
\param SolverHandle int Handle to an existing diffusion solver.

\return Success bool Returns if initialization failed.
\ingroup LuaDiffusion
\author Jan*/
int
chiDiffusionInitialize(lua_State* L)
{
  int solver_index = lua_tonumber(L, 1);

  auto& solver =
    Chi::GetStackItem<chi_diffusion::Solver>(Chi::object_stack, solver_index, __FUNCTION__);

  bool success = solver.Initialize(true);

  lua_pushnumber(L, success);
  return 1;
}
