#include "framework/lua.h"

#include "modules/DFEMDiffusion/dfem_diffusion_solver.h"

#include "framework/runtime.h"
#include "framework/logging/log.h"

int
chiDFEMDiffusionSolverCreate(lua_State* L)
{
  const std::string fname = __FUNCTION__;
  int num_args = lua_gettop(L);

  std::string solver_name = "DFEMDiffusionSolver";

  if (num_args == 1)
  {
    LuaCheckStringValue(fname, L, 1);
    solver_name = lua_tostring(L, 1);
  }

  auto new_solver = std::make_shared<dfem_diffusion::Solver>(solver_name);

  Chi::object_stack.push_back(new_solver);

  lua_pushinteger(L, static_cast<lua_Integer>(Chi::object_stack.size() - 1));

  Chi::log.LogAllVerbose1() << "\nchiDFEMDiffusionSolverCreate: DFEM Diffusion solver created"
                            << std::endl;
  return 1;
}
