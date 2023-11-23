#include "framework/lua.h"

#include "modules/diffusion_solver/diffusion_solver.h"

#include "framework/runtime.h"

#include "framework/runtime.h"
#include "framework/logging/log.h"

using namespace opensn;

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

  auto new_solver = std::make_shared<diffusion::Solver>(solver_name);

  opensn::object_stack.push_back(new_solver);

  lua_pushinteger(L, static_cast<lua_Integer>(opensn::object_stack.size() - 1));

  opensn::log.LogAllVerbose1() << "chiDiffusionCreateSolver: Diffusion solver created" << std::endl;
  return 1;
}
