#include "framework/chi_lua.h"

#include "modules/DiffusionSolver/Solver/diffusion_solver.h"

#include "framework/chi_runtime.h"

#include "framework/chi_runtime.h"
#include "framework/logging/chi_log.h"

/** Creates a Diffusion solver.

\return Handle int Handle to the created solver.
\ingroup LuaDiffusion
\author Jan*/
int
chiDiffusionCreateSolver(lua_State* L)
{
  const std::string fname = __FUNCTION__;
  int num_args = lua_gettop(L);

  std::string solver_name = "DiffusionSolver";

  if (num_args == 1)
  {
    LuaCheckStringValue(fname, L, 1);
    solver_name = lua_tostring(L, 1);
  }

  auto new_solver = std::make_shared<chi_diffusion::Solver>(solver_name);

  Chi::object_stack.push_back(new_solver);

  lua_pushinteger(L, static_cast<lua_Integer>(Chi::object_stack.size() - 1));

  Chi::log.LogAllVerbose1() << "chiDiffusionCreateSolver: Diffusion solver created" << std::endl;
  return 1;
}
