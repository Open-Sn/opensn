#include "framework/chi_lua.h"

#include "framework/mesh/MeshHandler/chi_meshhandler.h"
#include "framework/chi_runtime.h"

#include "framework/logging/chi_log.h"

#include <iostream>
#include "meshhandler_lua.h"
#include "framework/console/chi_console.h"

RegisterLuaFunctionAsIs(chiMeshHandlerCreate);
RegisterLuaFunctionAsIs(chiMeshHandlerSetCurrent);
RegisterLuaFunctionAsIs(chiMeshHandlerExportMeshToObj);
RegisterLuaFunctionAsIs(chiMeshHandlerExportMeshToVTK);
RegisterLuaFunctionAsIs(chiMeshHandlerExportMeshToExodus);

int
chiMeshHandlerCreate(lua_State* L)
{
  int index = (int)chi_mesh::PushNewHandlerAndGetIndex();
  lua_pushnumber(L, index);

  Chi::log.LogAllVerbose2() << "chiMeshHandlerCreate: Mesh Handler " << index << " created\n";

  return 1;
}

int
chiMeshHandlerSetCurrent(lua_State* L)
{
  int num_args = lua_gettop(L);
  if (num_args != 1) LuaPostArgAmountError("chiMeshHandlerSetCurrent", 1, num_args);

  int handle = lua_tonumber(L, 1);

  if ((handle < 0) or (handle >= Chi::meshhandler_stack.size()))
  {
    Chi::log.LogAllError() << "Invalid handle to mesh handler specified "
                           << "in call to chiMeshHandlerSetCurrent";
    Chi::Exit(EXIT_FAILURE);
  }

  Chi::current_mesh_handler = handle;

  Chi::log.LogAllVerbose2() << "chiMeshHandlerSetCurrent: set to " << handle;

  return 0;
}
