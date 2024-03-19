#include "lbs_common_lua_functions.h"
#include "modules/linear_boltzmann_solvers/lbs_solver/lbs_solver.h"
#include "framework/runtime.h"
#include "framework/lua.h"
#include "lua/framework/console/console.h"

using namespace opensn;

namespace opensnlua::lbs
{

RegisterLuaFunctionNamespace(LBSInitializeMaterials, lbs, InitializeMaterials);

int
LBSInitializeMaterials(lua_State* L)
{
  const std::string fname = "LBSInitializeMaterials";
  const int num_args = lua_gettop(L);

  if (num_args != 1)
    LuaPostArgAmountError(fname, 1, num_args);

  // Get pointer to solver
  const auto solver_handle = LuaArg<size_t>(L, 1);

  auto& lbs_solver =
    opensn::GetStackItem<opensn::lbs::LBSSolver>(opensn::object_stack, solver_handle, fname);

  lbs_solver.InitializeMaterials();

  return 0;
}

} // namespace opensnlua::lbs
