#include "framework/lua.h"

#include "modules/diffusion_solver/solver/diffusion_solver.h"

#include "framework/runtime.h"

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
