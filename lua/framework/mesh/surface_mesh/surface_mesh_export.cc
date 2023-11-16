#include "framework/lua.h"

#include <iostream>
#include "framework/mesh/surface_mesh/surface_mesh.h"
#include "framework/mesh/mesh_handler/mesh_handler.h"

#include "framework/runtime.h"

#include "framework/logging/log.h"
#include "lua_surface_mesh.h"
#include "framework/console/console.h"

using namespace opensn;

RegisterLuaFunctionAsIs(chiSurfaceMeshExportToObj);
RegisterLuaFunctionAsIs(chiSurfaceMeshExportPolyFile);

int
chiSurfaceMeshExportToObj(lua_State* L)
{
  auto& cur_hndlr = opensn::GetCurrentHandler();

  // Get arguments
  int num_args = lua_gettop(L);
  if (num_args != 2) LuaPostArgAmountError("chiSurfaceMeshExportObj", 2, num_args);

  int handle = lua_tonumber(L, 1);

  size_t length = 0;
  const char* temp = lua_tolstring(L, 2, &length);

  auto& surface_mesh =
    opensn::Chi::GetStackItem<SurfaceMesh>(opensn::Chi::surface_mesh_stack, handle, __FUNCTION__);

  surface_mesh.ExportToOBJFile(temp);

  return 0;
}

int
chiSurfaceMeshExportPolyFile(lua_State* L)
{
  auto& cur_hndlr = opensn::GetCurrentHandler();

  // Get arguments
  int num_args = lua_gettop(L);
  if (num_args != 2) LuaPostArgAmountError("chiSurfaceMeshExportPolyFile", 2, num_args);

  int handle = lua_tonumber(L, 1);

  size_t length = 0;
  const char* temp = lua_tolstring(L, 2, &length);

  auto& surface_mesh =
    opensn::Chi::GetStackItem<SurfaceMesh>(opensn::Chi::surface_mesh_stack, handle, __FUNCTION__);

  surface_mesh.ExportToPolyFile(temp);
  return 0;
}
