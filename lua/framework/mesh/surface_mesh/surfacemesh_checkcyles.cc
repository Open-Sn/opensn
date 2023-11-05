#include "opensn/framework/chi_lua.h"

#include <iostream>
#include <algorithm>
#include "opensn/framework/mesh/SurfaceMesh/chi_surfacemesh.h"
#include "opensn/framework/mesh/MeshHandler/chi_meshhandler.h"
#include "opensn/framework/chi_runtime.h"

#include "opensn/framework/logging/chi_log.h"
#include "lua_surface_mesh.h"
#include "opensn/framework/console/chi_console.h"

RegisterLuaFunctionAsIs(chiSurfaceMeshCheckCycles);
RegisterLuaFunctionAsIs(chiComputeLoadBalancing);

//#############################################################################
/** Builds sweep ordering for a number of angles and checks whether any
 * cyclic dependencies are encountered.
 *
\param SurfaceHandle int Handle to the surface on which the operation is to be
performed. \param NumAngles int Number of azimuthal angles to use for checking
cycles.

\ingroup LuaSurfaceMesh
\author Jan*/
int
chiSurfaceMeshCheckCycles(lua_State* L)
{
  int num_args = lua_gettop(L);
  if (num_args != 2) LuaPostArgAmountError("chiSurfaceMeshCheckCycles", 2, num_args);

  auto& cur_hndlr = chi_mesh::GetCurrentHandler();

  int surf_handle = lua_tonumber(L, 1);
  int num_angles = lua_tonumber(L, 2);

  auto& surf_mesh =
    Chi::GetStackItem<chi_mesh::SurfaceMesh>(Chi::surface_mesh_stack, surf_handle, __FUNCTION__);

  surf_mesh.CheckCyclicDependencies(num_angles);
  return 0;
}

//#############################################################################
/** Computes load balancing parameters for given predictive x and y cuts
 * without actually performing cuts.
 *
\param SurfaceHandle int Handle to the surface on which the operation is to be
performed. \param Xcuts table Array of x-values associated with the xcuts.
\param Ycuts table Array of y-values associated with the ycuts.

\ingroup LuaSurfaceMesh
\author Jan*/
int
chiComputeLoadBalancing(lua_State* L)
{
  int num_args = lua_gettop(L);
  if (num_args != 3) LuaPostArgAmountError("chiComputeLoadBalancing", 3, num_args);

  //======================================== Get reference surface mesh
  int surf_handle = lua_tonumber(L, 1);

  auto& cur_surf =
    Chi::GetStackItem<chi_mesh::SurfaceMesh>(Chi::surface_mesh_stack, surf_handle, __FUNCTION__);

  //======================================== Extract x-cuts
  if (!lua_istable(L, 2))
  {
    Chi::log.LogAllError() << "In call to chiComputeLoadBalancing: "
                           << " expected table for argument 2. Incompatible value supplied.";
    Chi::Exit(EXIT_FAILURE);
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

  //======================================== Extract y-cuts
  if (!lua_istable(L, 3))
  {
    Chi::log.LogAllError() << "In call to chiComputeLoadBalancing: "
                           << " expected table for argument 3. Incompatible value supplied.";
    Chi::Exit(EXIT_FAILURE);
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

  //======================================== Call compute balance
  std::stable_sort(x_cuts.begin(), x_cuts.end());
  std::stable_sort(y_cuts.begin(), y_cuts.end());
  cur_surf.ComputeLoadBalancing(x_cuts, y_cuts);

  return 0;
}
