#include "modules/linear_boltzmann_solvers/a_lbs_solver/lbs_solver.h"

#include "framework/runtime.h"
#include "framework/logging/log.h"
#include "modules/linear_boltzmann_solvers/a_lbs_solver/groupset/lbs_groupset.h"
#include "framework/lua.h"

using namespace opensn;

namespace opensnlua::lbs
{

int
chiLBSWriteGroupsetAngularFlux(lua_State* L)
{
  const std::string fname = "chiLBSWriteGroupsetAngularFlux";
  // Get arguments
  const int num_args = lua_gettop(L);
  if (num_args != 3) LuaPostArgAmountError(fname, 3, num_args);

  LuaCheckNilValue(fname, L, 1);
  LuaCheckNilValue(fname, L, 2);
  LuaCheckNilValue(fname, L, 3);

  const int solver_handle = lua_tonumber(L, 1);
  const int grpset_index = lua_tonumber(L, 2);
  const std::string file_base = lua_tostring(L, 3);

  // Get pointer to solver
  auto& lbs_solver =
    opensn::GetStackItem<opensn::lbs::LBSSolver>(opensn::object_stack, solver_handle, fname);

  // Obtain pointer to groupset
  opensn::lbs::LBSGroupset* groupset = nullptr;
  try
  {
    groupset = &lbs_solver.Groupsets().at(grpset_index);
  }
  catch (const std::out_of_range& o)
  {
    opensn::log.LogAllError() << "Invalid handle to groupset "
                              << "in call to " << fname;
    opensn::Exit(EXIT_FAILURE);
  }

  const auto& psi = lbs_solver.PsiNewLocal().at(groupset->id_);
  lbs_solver.WriteGroupsetAngularFluxes(*groupset, psi, file_base);

  return 0;
}

int
chiLBSReadGroupsetAngularFlux(lua_State* L)
{
  const std::string fname = "chiLBSReadGroupsetAngularFlux";
  // Get arguments
  const int num_args = lua_gettop(L);
  if (num_args != 3) LuaPostArgAmountError(fname, 3, num_args);

  LuaCheckNilValue(fname, L, 1);
  LuaCheckNilValue(fname, L, 2);
  LuaCheckNilValue(fname, L, 3);

  const int solver_handle = lua_tonumber(L, 1);
  const int grpset_index = lua_tonumber(L, 2);
  const std::string file_base = lua_tostring(L, 3);

  // Get pointer to solver
  auto& lbs_solver =
    opensn::GetStackItem<opensn::lbs::LBSSolver>(opensn::object_stack, solver_handle, fname);

  // Obtain pointer to groupset
  opensn::lbs::LBSGroupset* groupset = nullptr;
  try
  {
    groupset = &lbs_solver.Groupsets().at(grpset_index);
  }
  catch (const std::out_of_range& o)
  {
    opensn::log.LogAllError() << "Invalid handle to groupset "
                              << "in call to " << fname;
    opensn::Exit(EXIT_FAILURE);
  }

  auto& psi = lbs_solver.PsiNewLocal().at(groupset->id_);
  lbs_solver.ReadGroupsetAngularFluxes(file_base, *groupset, psi);

  return 0;
}

} // namespace opensnlua::lbs
