// SPDX-FileCopyrightText: 2024 The OpenSn Authors <https://open-sn.github.io/opensn/>
// SPDX-License-Identifier: MIT

#include "lua/modules/linear_bolzmann_solvers/lbs_solver/lbs_lua_utils.h"
#include "lua/framework/math/functions/lua_vector_spatial_material_function.h"
#include "lua/framework/console/console.h"
#include "lua/framework/lua.h"
#include "framework/runtime.h"
#include "framework/logging/log.h"
#include "framework/mesh/logical_volume/logical_volume.h"
#include "modules/linear_boltzmann_solvers/lbs_solver/lbs_solver.h"

using namespace opensn;

namespace opensnlua
{

RegisterLuaFunctionInNamespace(AddPointSource, lbs, AddPointSource);
int
AddPointSource(lua_State* L)
{
  const std::string fname = "lbs.AddPointSource";
  LuaCheckArgs<int, int>(L, fname);

  // Process solver handle
  const auto solver_handle = LuaArg<int>(L, 1);
  auto& lbs_solver =
    opensn::GetStackItem<opensn::LBSSolver>(opensn::object_stack, solver_handle, __FUNCTION__);

  // Process point source handle
  auto pt_src_handle = LuaArg<int>(L, 2);
  lbs_solver.AddPointSource(std::move(
    opensn::GetStackItem<opensn::PointSource>(opensn::object_stack, pt_src_handle, __FUNCTION__)));
  return LuaReturn(L);
}

RegisterLuaFunctionInNamespace(ClearPointSources, lbs, ClearPointSources);
int
ClearPointSources(lua_State* L)
{
  const std::string fname = "lbs.ClearPointSources";
  LuaCheckArgs<int>(L, fname);

  const auto solver_handle = LuaArg<int>(L, 1);
  auto& lbs_solver =
    opensn::GetStackItem<opensn::LBSSolver>(opensn::object_stack, solver_handle, __FUNCTION__);

  lbs_solver.ClearPointSources();
  opensn::log.Log() << "Cleared all point sources.";
  return LuaReturn(L);
}

RegisterLuaFunctionInNamespace(AddVolumetricSource, lbs, AddVolumetricSource);
int
AddVolumetricSource(lua_State* L)
{
  const std::string fname = "lbs.AddVolumetricSource";
  LuaCheckArgs<int, int>(L, fname);

  const auto solver_handle = LuaArg<int>(L, 1);
  auto& lbs_solver =
    opensn::GetStackItem<opensn::LBSSolver>(object_stack, solver_handle, __FUNCTION__);

  const auto src_handle = LuaArg<int>(L, 2);
  lbs_solver.AddVolumetricSource(std::move(
    opensn::GetStackItem<opensn::VolumetricSource>(object_stack, src_handle, __FUNCTION__)));

  opensn::log.Log() << lbs_solver.Name() << ": Added volumetric source.";
  return LuaReturn(L);
}

RegisterLuaFunctionInNamespace(ClearVolumetricSources, lbs, ClearVolumetricSources);
int
ClearVolumetricSources(lua_State* L)
{
  const std::string fname = "lbs.ClearVolumetricSources";
  LuaCheckArgs<int>(L, fname);

  // Process solver handle
  const auto solver_handle = LuaArg<int>(L, 1);
  auto& lbs_solver =
    opensn::GetStackItem<opensn::LBSSolver>(opensn::object_stack, solver_handle, __FUNCTION__);

  lbs_solver.ClearVolumetricSources();
  opensn::log.Log() << "Cleared all volumetric sources.";
  return LuaReturn(L);
}

} // namespace opensnlua
