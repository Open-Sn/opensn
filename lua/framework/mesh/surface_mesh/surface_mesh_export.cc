#include "framework/lua.h"

#include <iostream>
#include "framework/mesh/surface_mesh/surface_mesh.h"

#include "framework/runtime.h"

#include "framework/logging/log.h"
#include "lua_surface_mesh.h"
#include "framework/console/console.h"

using namespace opensn;

RegisterLuaFunctionAsIs(SurfaceMeshExportToObj);
RegisterLuaFunctionAsIs(SurfaceMeshExportPolyFile);

int
SurfaceMeshExportToObj(lua_State* L)
{
  // Get arguments
  int num_args = lua_gettop(L);
  if (num_args != 2)
    LuaPostArgAmountError("SurfaceMeshExportObj", 2, num_args);

  int handle = lua_tonumber(L, 1);

  size_t length = 0;
  const char* temp = lua_tolstring(L, 2, &length);

  auto& surface_mesh =
    opensn::GetStackItem<SurfaceMesh>(opensn::surface_mesh_stack, handle, __FUNCTION__);

  surface_mesh.ExportToOBJFile(temp);

  return 0;
}

int
SurfaceMeshExportPolyFile(lua_State* L)
{
  // Get arguments
  int num_args = lua_gettop(L);
  if (num_args != 2)
    LuaPostArgAmountError("SurfaceMeshExportPolyFile", 2, num_args);

  int handle = lua_tonumber(L, 1);

  size_t length = 0;
  const char* temp = lua_tolstring(L, 2, &length);

  auto& surface_mesh =
    opensn::GetStackItem<SurfaceMesh>(opensn::surface_mesh_stack, handle, __FUNCTION__);

  surface_mesh.ExportToPolyFile(temp);
  return 0;
}
