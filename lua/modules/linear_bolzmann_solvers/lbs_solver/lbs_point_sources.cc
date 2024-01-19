#include "lbs_lua_utils.h"

#include "framework/runtime.h"
#include "framework/console/console.h"
#include "framework/logging/log.h"
#include "framework/lua.h"

#include "modules/linear_boltzmann_solvers/a_lbs_solver/lbs_solver.h"

using namespace opensn;

namespace opensnlua::lbs
{

int
chiLBSAddPointSource(lua_State* L)
{
  opensn::log.Log0Warning() << "chiLBSAddPointSource has been deprecated and will be removed soon. "
                               "Consider using lbs.AddPointSource instead, or setting point "
                               "sources via LBSSolver::Options.";

  const std::string fname = __FUNCTION__;

  const int num_args = lua_gettop(L);
  if (num_args != 5) LuaPostArgAmountError(fname, 5, num_args);

  LuaCheckNilValue(fname, L, 1);
  LuaCheckNilValue(fname, L, 2);
  LuaCheckNilValue(fname, L, 3);
  LuaCheckNilValue(fname, L, 4);
  LuaCheckNilValue(fname, L, 5);

  // Process handle
  LuaCheckIntegerValue(fname, L, 1);
  const int solver_handle = lua_tointeger(L, 1);
  auto& lbs_solver = opensn::Chi::GetStackItem<opensn::lbs::LBSSolver>(
    opensn::Chi::object_stack, solver_handle, fname);

  // Proces location
  LuaCheckNumberValue(fname, L, 2);
  LuaCheckNumberValue(fname, L, 3);
  LuaCheckNumberValue(fname, L, 4);

  const double x = lua_tonumber(L, 2);
  const double y = lua_tonumber(L, 3);
  const double z = lua_tonumber(L, 4);
  const std::vector<double> location = {x, y, z};

  // Process source strength
  LuaCheckTableValue(fname, L, 5);
  std::vector<double> strength;
  LuaPopulateVectorFrom1DArray(fname, L, 5, strength);

  // Create point source specification
  ParameterBlock params;
  params.AddParameter("location", location);
  params.AddParameter("strength", strength);

  auto pt_src_spec = opensn::lbs::PointSource::GetInputParameters();
  pt_src_spec.AssignParameters(params);

  lbs_solver.AddPointSource(opensn::lbs::PointSource(pt_src_spec));

  opensn::log.Log() << "Added point source at location "
                    << lbs_solver.PointSources().back().Location().PrintStr();

  return 0;
}

int
chiLBSClearPointSources(lua_State* L)
{
  opensn::log.Log0Warning() << "chiLBSClearPointSource has been deprecated and will "
                               "be removed soon. Consider using lbs.ClearPointSources instead, "
                               "or clearing point sources via LBSSolver::Options.";

  const std::string fname = __FUNCTION__;

  const int num_args = lua_gettop(L);
  if (num_args != 1) LuaPostArgAmountError(fname, 1, num_args);

  // Process handle
  LuaCheckIntegerValue(fname, L, 1);
  const int solver_handle = lua_tointeger(L, 1);
  auto& lbs_solver = opensn::Chi::GetStackItem<opensn::lbs::LBSSolver>(
    opensn::Chi::object_stack, solver_handle, fname);

  lbs_solver.ClearPointSources();

  opensn::log.Log() << "Cleared all point sources.";

  return 0;
}

RegisterLuaFunction(AddPointSource, lbs, AddPointSource);
int
AddPointSource(lua_State* L)
{
  const std::string fname = __FUNCTION__;

  const int num_args = lua_gettop(L);
  if (num_args != 2) LuaPostArgAmountError(fname, 2, num_args);

  // Process handle
  LuaCheckIntegerValue(fname, L, 1);
  const int solver_handle = lua_tointeger(L, 1);
  auto& lbs_solver = opensn::Chi::GetStackItem<opensn::lbs::LBSSolver>(
    opensn::Chi::object_stack, solver_handle, fname);

  // Process point source specification
  LuaCheckIntegerValue(fname, L, 2);
  int pt_src_handle = lua_tointeger(L, 2);
  lbs_solver.AddPointSource(opensn::Chi::GetStackItem<opensn::lbs::PointSource>(
    opensn::Chi::object_stack, pt_src_handle, fname));
  return 1;
}

RegisterLuaFunction(ClearPointSources, lbs, ClearPointSources);
int
ClearPointSources(lua_State* L)
{
  const std::string fname = __FUNCTION__;

  const int num_args = lua_gettop(L);
  if (num_args != 1) LuaPostArgAmountError(fname, 1, num_args);

  // Process handle
  LuaCheckIntegerValue(fname, L, 1);
  const int solver_handle = lua_tointeger(L, 1);
  auto& lbs_solver = opensn::Chi::GetStackItem<opensn::lbs::LBSSolver>(
    opensn::Chi::object_stack, solver_handle, fname);

  lbs_solver.ClearPointSources();
  opensn::log.Log() << "Cleared all point sources.";
  return 0;
}

} // namespace opensnlua::lbs
