#include "lbsadj_lua_utils.h"

#include "modules/linear_boltzmann_solvers/discrete_ordinates_adjoint_solver/lbs_adj_solver.h"

#include "framework/runtime.h"

#include "framework/console/console.h"

using namespace opensn;

namespace opensnlua::lbs
{

RegisterLuaFunctionAsIs(AdjointSolverExportImportanceMapBinary);

int
AdjointSolverExportImportanceMapBinary(lua_State* L)
{
  const std::string fname = __FUNCTION__;
  const int num_args = lua_gettop(L);
  if (num_args != 2)
    LuaPostArgAmountError(fname, 2, num_args);

  const auto solver_handle = LuaArg<int>(L, 1);
  const auto file_name = LuaArg<std::string>(L, 2);

  auto& solver = opensn::GetStackItem<opensn::lbs::DiscreteOrdinatesAdjointSolver>(
    opensn::object_stack, solver_handle, fname);

  solver.ExportImportanceMap(file_name);

  return 0;
}

} // namespace opensnlua::lbs
