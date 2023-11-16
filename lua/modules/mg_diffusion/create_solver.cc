#include "framework/lua.h"

#include "modules/mg_diffusion/mg_diffusion_solver.h"

#include "framework/runtime.h"
#include "framework/logging/log.h"

using namespace opensn;

namespace opensnlua::mg_diffusion
{

int
chiCFEMMGDiffusionSolverCreate(lua_State* L)
{
  const std::string fname = __FUNCTION__;
  int num_args = lua_gettop(L);

  std::string solver_name = "MGDiffusionSolver";

  if (num_args == 1)
  {
    LuaCheckStringValue(fname, L, 1);
    solver_name = lua_tostring(L, 1);
  }

  auto new_solver = std::make_shared<opensn::mg_diffusion::Solver>(solver_name);

  opensn::Chi::object_stack.push_back(new_solver);

  lua_pushinteger(L, static_cast<lua_Integer>(opensn::Chi::object_stack.size() - 1));

  opensn::Chi::log.LogAllVerbose1() << "\nchiCFEMMGDiffusionSolverCreate: CFEM "
                                       "Multigroup Diffusion solver created"
                                    << std::endl;
  return 1;
}

} // namespace opensnlua::mg_diffusion
