#include "modules/linear_boltzmann_solvers/a_lbs_solver/lbs_solver.h"

#include "framework/runtime.h"
#include "framework/lua.h"

using namespace opensn;

namespace opensnlua::lbs
{

int
LBSInitializeMaterials(lua_State* L)
{
  const std::string fname = "LBSInitializeMaterials";
  const int num_args = lua_gettop(L);

  if (num_args != 1)
    LuaPostArgAmountError(fname, 1, num_args);

  LuaCheckNilValue(fname, L, 1);

  // Get pointer to solver
  const int solver_handle = lua_tonumber(L, 1);

  auto& lbs_solver =
    opensn::GetStackItem<opensn::lbs::LBSSolver>(opensn::object_stack, solver_handle, fname);

  lbs_solver.InitializeMaterials();

  return 0;
}

} // namespace opensnlua::lbs
