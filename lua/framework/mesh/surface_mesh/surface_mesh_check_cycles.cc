#include "lua_surface_mesh.h"
#include "framework/mesh/surface_mesh/surface_mesh.h"
#include "framework/runtime.h"
#include "framework/logging/log.h"
#include "framework/console/console.h"
#include <iostream>
#include <algorithm>

using namespace opensn;

RegisterLuaFunctionNamespace(MeshComputeLoadBalancing, mesh, ComputeLoadBalancing);

int
MeshComputeLoadBalancing(lua_State* L)
{
  int num_args = lua_gettop(L);
  if (num_args != 3)
    LuaPostArgAmountError("mesh.ComputeLoadBalancing", 3, num_args);

  // Get reference surface mesh
  int surf_handle = lua_tonumber(L, 1);

  auto& cur_surf =
    opensn::GetStackItem<SurfaceMesh>(opensn::surface_mesh_stack, surf_handle, __FUNCTION__);

  // Extract x-cuts
  if (not lua_istable(L, 2))
  {
    opensn::log.LogAllError() << "In call to mesh.ComputeLoadBalancing: expected table for "
                                 "argument 2. Incompatible value supplied.";
    opensn::Exit(EXIT_FAILURE);
  }

  int x_table_len = lua_rawlen(L, 2);

  std::vector<double> x_cuts(x_table_len, 0.0);
  for (int g = 0; g < x_table_len; g++)
  {
    lua_pushnumber(L, g + 1);
    lua_gettable(L, 2);
    x_cuts[g] = lua_tonumber(L, -1);
    lua_pop(L, 1);
  }

  // Extract y-cuts
  if (not lua_istable(L, 3))
  {
    opensn::log.LogAllError() << "In call to mesh.ComputeLoadBalancing: expected table for "
                                 "argument 3. Incompatible value supplied.";
    opensn::Exit(EXIT_FAILURE);
  }

  int y_table_len = lua_rawlen(L, 3);

  std::vector<double> y_cuts(y_table_len, 0.0);
  for (int g = 0; g < y_table_len; g++)
  {
    lua_pushnumber(L, g + 1);
    lua_gettable(L, 3);
    y_cuts[g] = lua_tonumber(L, -1);
    lua_pop(L, 1);
  }

  // Call compute balance
  std::stable_sort(x_cuts.begin(), x_cuts.end());
  std::stable_sort(y_cuts.begin(), y_cuts.end());
  cur_surf.ComputeLoadBalancing(x_cuts, y_cuts);

  return 0;
}
