#include "lbs_do_lua_utils.h"

#include "modules/linear_boltzmann_solvers/b_discrete_ordinates_solver/lbs_discrete_ordinates_solver.h"

#include "framework/console/console.h"

#include "framework/runtime.h"

using namespace opensn;

namespace opensnlua::lbs
{

RegisterLuaFunctionAsIs(chiLBSComputeLeakage);

int
chiLBSComputeLeakage(lua_State* L)
{
  const std::string fname = "chiLBSComputeLeakage";
  const int num_args = lua_gettop(L);

  if (num_args != 3) LuaPostArgAmountError(fname, 3, num_args);

  LuaCheckNilValue(fname, L, 1);

  // Get pointer to solver
  const int solver_handle = lua_tonumber(L, 1);

  auto& lbs_solver = opensn::Chi::GetStackItem<opensn::lbs::DiscreteOrdinatesSolver>(
    opensn::Chi::object_stack, solver_handle, fname);
  LuaCheckNilValue(fname, L, 2);
  LuaCheckNilValue(fname, L, 3);

  const int groupset_id = lua_tonumber(L, 2);
  const int boundary_id = lua_tonumber(L, 3);

  const auto leakage = lbs_solver.ComputeLeakage(groupset_id, boundary_id);

  // Push up the table
  lua_newtable(L);

  for (int i = 0; i < static_cast<int>(leakage.size()); ++i)
  {
    lua_pushinteger(L, i + 1);
    lua_pushnumber(L, leakage[i]);
    lua_settable(L, -3);
  }

  return 1;
}

} // namespace opensnlua::lbs
