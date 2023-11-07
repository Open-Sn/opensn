#include "modules/LinearBoltzmannSolvers/A_LBSSolver/lbs_solver.h"

#include "framework/runtime.h"

namespace lbs::common_lua_utils
{

int
chiLBSInitializeMaterials(lua_State* L)
{
  const std::string fname = "chiLBSInitializeMaterials";
  const int num_args = lua_gettop(L);

  if (num_args != 1) LuaPostArgAmountError(fname, 1, num_args);

  LuaCheckNilValue(fname, L, 1);

  // Get pointer to solver
  const int solver_handle = lua_tonumber(L, 1);

  auto& lbs_solver = Chi::GetStackItem<lbs::LBSSolver>(Chi::object_stack, solver_handle, fname);

  lbs_solver.InitMaterials();

  return 0;
}

} // namespace lbs::common_lua_utils
