#include "framework/lua.h"
#include "modules/diffusion_solver/diffusion_solver.h"

#include "framework/runtime.h"

using namespace opensn;

int
DiffusionExecute(lua_State* L)
{
  const size_t solver_index = lua_tonumber(L, 1);
  auto& solver =
    opensn::GetStackItem<diffusion::Solver>(opensn::object_stack, solver_index, __FUNCTION__);

  solver.ExecuteS();

  return 0;
}
