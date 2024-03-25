#include "lbs_common_lua_functions.h"
#include "modules/linear_boltzmann_solvers/lbs_solver/lbs_solver.h"
#include "framework/runtime.h"
#include "framework/logging/log.h"
#include "modules/linear_boltzmann_solvers/lbs_solver/groupset/lbs_groupset.h"
#include "framework/lua.h"
#include "lua/framework/console/console.h"

using namespace opensn;

namespace opensnlua::lbs
{

RegisterLuaFunctionNamespace(LBSWriteGroupsetAngularFlux, lbs, WriteGroupsetAngularFlux);
RegisterLuaFunctionNamespace(LBSReadGroupsetAngularFlux, lbs, ReadGroupsetAngularFlux);

int
LBSWriteGroupsetAngularFlux(lua_State* L)
{
  const std::string fname = "lbs.WriteGroupsetAngularFlux";
  LuaCheckArgs<size_t, int, std::string>(L, fname);

  const auto solver_handle = LuaArg<size_t>(L, 1);
  const auto grpset_index = LuaArg<int>(L, 2);
  const auto file_base = LuaArg<std::string>(L, 3);

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

  return LuaReturn(L);
}

int
LBSReadGroupsetAngularFlux(lua_State* L)
{
  const std::string fname = "lbs.ReadGroupsetAngularFlux";
  LuaCheckArgs<size_t, int, std::string>(L, fname);

  const auto solver_handle = LuaArg<size_t>(L, 1);
  const auto grpset_index = LuaArg<int>(L, 2);
  const auto file_base = LuaArg<std::string>(L, 3);

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

  return LuaReturn(L);
}

} // namespace opensnlua::lbs
