#include "framework/lua.h"
#include "modules/diffusion_solver/diffusion_solver.h"

#include "framework/runtime.h"

using namespace opensn;

int
chiDiffusionExecute(lua_State* L)
{
  const size_t solver_index = lua_tonumber(L, 1);
  auto& solver = opensn::Chi::GetStackItem<diffusion::Solver>(
    opensn::Chi::object_stack, solver_index, __FUNCTION__);

  solver.ExecuteS();

  return 0;
}
