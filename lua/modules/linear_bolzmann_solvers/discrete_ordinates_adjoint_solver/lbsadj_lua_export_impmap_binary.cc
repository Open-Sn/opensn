#include "lbsadj_lua_utils.h"

#include "modules/LinearBoltzmannSolvers/C_DiscreteOrdinatesAdjointSolver/lbsadj_solver.h"

#include "framework/runtime.h"

#include "framework/console/console.h"

namespace lbs::adjoint_lua_utils
{

RegisterLuaFunctionAsIs(chiAdjointSolverExportImportanceMapBinary);

int
chiAdjointSolverExportImportanceMapBinary(lua_State* L)
{
  const std::string fname = __FUNCTION__;
  const int num_args = lua_gettop(L);
  if (num_args != 2) LuaPostArgAmountError(fname, 2, num_args);

  LuaCheckNilValue(fname, L, 1);
  LuaCheckNilValue(fname, L, 2);

  LuaCheckIntegerValue(fname, L, 1);
  LuaCheckStringValue(fname, L, 2);

  const int solver_handle = lua_tointeger(L, 1);
  const std::string file_name = lua_tostring(L, 2);

  auto& solver =
    Chi::GetStackItem<lbs::DiscreteOrdinatesAdjointSolver>(Chi::object_stack, solver_handle, fname);

  solver.ExportImportanceMap(file_name);

  return 0;
}

} // namespace lbs::adjoint_lua_utils
