#include "chi_lua.h"

#include"../mg_diffusion_solver.h"

#include "chi_runtime.h"
#include "chi_log.h"

namespace mg_diffusion::mgd_lua_utils
{

//#############################################################################
/** Creates a Multigroup CFEM Diffusion solver.

\return Handle int Handle to the created solver.
\ingroup LuaDiffusion
*/
int chiCFEMMGDiffusionSolverCreate(lua_State *L)
{
  const std::string fname = __FUNCTION__;
  int num_args = lua_gettop(L);

  std::string solver_name = "MGDiffusionSolver";

  if (num_args == 1)
  {
    LuaCheckStringValue(fname, L, 1);
    solver_name = lua_tostring(L, 1);
  }

  auto new_solver = std::make_shared<mg_diffusion::Solver>(solver_name);

  Chi::object_stack.push_back(new_solver);

  lua_pushinteger(L,
      static_cast<lua_Integer>(Chi::object_stack.size()-1));

  Chi::log.LogAllVerbose1()
    << "\nchiCFEMMGDiffusionSolverCreate: CFEM Multigroup Diffusion solver created"
    << std::endl;
  return 1;
}

}//namespace mg_diffusion::mgd_lua_utils