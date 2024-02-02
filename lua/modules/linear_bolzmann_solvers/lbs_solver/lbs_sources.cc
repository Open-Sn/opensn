#include "lbs_lua_utils.h"

#include "framework/runtime.h"
#include "framework/console/console.h"
#include "framework/logging/log.h"
#include "framework/lua.h"

#include "framework/mesh/logical_volume/logical_volume.h"
#include "lua/framework/math/functions/lua_spatial_material_function.h"
#include "modules/linear_boltzmann_solvers/a_lbs_solver/lbs_solver.h"

using namespace opensn;

namespace opensnlua::lbs
{

int
LBSAddPointSource(lua_State* L)
{
  opensn::log.Log0Warning() << "LBSAddPointSource has been deprecated and will be removed soon. "
                               "Consider using lbs.AddPointSource instead, or setting point "
                               "sources via LBSSolver::Options.";

  const int num_args = lua_gettop(L);
  if (num_args != 5) LuaPostArgAmountError(__FUNCTION__, 5, num_args);

  LuaCheckIntegerValue(__FUNCTION__, L, 1);
  LuaCheckNumberValue(__FUNCTION__, L, 2);
  LuaCheckNumberValue(__FUNCTION__, L, 3);
  LuaCheckNumberValue(__FUNCTION__, L, 4);
  LuaCheckTableValue(__FUNCTION__, L, 5);

  // Process solver handle
  const int solver_handle = lua_tointeger(L, 1);
  auto& lbs_solver =
    opensn::GetStackItem<opensn::lbs::LBSSolver>(opensn::object_stack, solver_handle, __FUNCTION__);

  // Proces location
  const double x = lua_tonumber(L, 2);
  const double y = lua_tonumber(L, 3);
  const double z = lua_tonumber(L, 4);
  const std::vector<double> location = {x, y, z};

  // Process source strength
  std::vector<double> strength;
  LuaPopulateVectorFrom1DArray(__FUNCTION__, L, 5, strength);

  // Create point source specification
  ParameterBlock params;
  params.AddParameter("location", location);
  params.AddParameter("strength", strength);

  auto spec = opensn::lbs::PointSource::GetInputParameters();
  spec.AssignParameters(params);

  lbs_solver.AddPointSource(opensn::lbs::PointSource(spec));
  opensn::log.Log() << "Added point source at location "
                    << lbs_solver.PointSources().back().Location().PrintStr();
  return 0;
}

int
LBSClearPointSources(lua_State* L)
{
  opensn::log.Log0Warning() << "chiLBSClearPointSource has been deprecated and will "
                               "be removed soon. Consider using lbs.ClearPointSources instead, "
                               "or clearing point sources via LBSSolver::Options.";

  const int num_args = lua_gettop(L);
  if (num_args != 1) LuaPostArgAmountError(__FUNCTION__, 1, num_args);

  LuaCheckIntegerValue(__FUNCTION__, L, 1);

  // Process solver handle
  const int solver_handle = lua_tointeger(L, 1);
  auto& lbs_solver =
    opensn::GetStackItem<opensn::lbs::LBSSolver>(opensn::object_stack, solver_handle, __FUNCTION__);

  lbs_solver.ClearPointSources();
  opensn::log.Log() << "Cleared all point sources.";
  return 0;
}

RegisterLuaFunctionNamespace(AddPointSource, lbs, AddPointSource);
int
AddPointSource(lua_State* L)
{
  const int num_args = lua_gettop(L);
  if (num_args != 2) LuaPostArgAmountError(__FUNCTION__, 2, num_args);

  LuaCheckIntegerValue(__FUNCTION__, L, 1);
  LuaCheckIntegerValue(__FUNCTION__, L, 2);

  // Process solver handle
  const int solver_handle = lua_tointeger(L, 1);
  auto& lbs_solver =
    opensn::GetStackItem<opensn::lbs::LBSSolver>(opensn::object_stack, solver_handle, __FUNCTION__);

  // Process point source handle
  int pt_src_handle = lua_tointeger(L, 2);
  lbs_solver.AddPointSource(std::move(opensn::GetStackItem<opensn::lbs::PointSource>(
    opensn::object_stack, pt_src_handle, __FUNCTION__)));
  return 1;
}

RegisterLuaFunctionNamespace(ClearPointSources, lbs, ClearPointSources);
int
ClearPointSources(lua_State* L)
{
  const int num_args = lua_gettop(L);
  if (num_args != 1) LuaPostArgAmountError(__FUNCTION__, 1, num_args);

  LuaCheckIntegerValue(__FUNCTION__, L, 1);

  // Process solver handle
  const int solver_handle = lua_tointeger(L, 1);
  auto& lbs_solver =
    opensn::GetStackItem<opensn::lbs::LBSSolver>(opensn::object_stack, solver_handle, __FUNCTION__);

  lbs_solver.ClearPointSources();
  opensn::log.Log() << "Cleared all point sources.";
  return 0;
}

RegisterLuaFunctionNamespace(AddDistributedSource, lbs, AddDistributedSource);
int
AddDistributedSource(lua_State* L)
{
  const int num_args = lua_gettop(L);
  if (num_args != 2) LuaPostArgAmountError(__FUNCTION__, 2, num_args);

  LuaCheckIntegerValue(__FUNCTION__, L, 1);
  LuaCheckIntegerValue(__FUNCTION__, L, 2);

  const int solver_handle = lua_tointeger(L, 1);
  auto& lbs_solver =
    opensn::GetStackItem<opensn::lbs::LBSSolver>(object_stack, solver_handle, __FUNCTION__);

  const int src_handle = lua_tointeger(L, 2);
  lbs_solver.AddDistributedSource(std::move(
    opensn::GetStackItem<opensn::lbs::DistributedSource>(object_stack, src_handle, __FUNCTION__)));

  opensn::log.Log() << lbs_solver.TextName() << ": Added distributed source.";
  return 0;
}

RegisterLuaFunctionNamespace(ClearDistributedSources, lbs, ClearDistributedSources);
int
ClearDistributedSources(lua_State* L)
{
  const int num_args = lua_gettop(L);
  if (num_args != 1) LuaPostArgAmountError(__FUNCTION__, 1, num_args);

  LuaCheckIntegerValue(__FUNCTION__, L, 1);

  // Process solver handle
  const int solver_handle = lua_tointeger(L, 1);
  auto& lbs_solver =
    opensn::GetStackItem<opensn::lbs::LBSSolver>(opensn::object_stack, solver_handle, __FUNCTION__);

  lbs_solver.ClearDistributedSources();
  opensn::log.Log() << "Cleared all distributed sources.";
  return 0;
}

} // namespace opensnlua::lbs
