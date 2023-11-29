#include "framework/lua.h"

#include "framework/mesh/mesh_handler/mesh_handler.h"
#include "framework/runtime.h"

#include "framework/logging/log.h"

#include <iostream>
#include "mesh_handler_lua.h"
#include "framework/console/console.h"

using namespace opensn;

RegisterLuaFunctionAsIs(MeshHandlerCreate);
RegisterLuaFunctionAsIs(chiMeshHandlerSetCurrent);
RegisterLuaFunctionAsIs(MeshHandlerExportMeshToObj);
RegisterLuaFunctionAsIs(MeshHandlerExportMeshToVTK);
RegisterLuaFunctionAsIs(MeshHandlerExportMeshToExodus);

int
MeshHandlerCreate(lua_State* L)
{
  int index = (int)PushNewHandlerAndGetIndex();
  lua_pushnumber(L, index);

  opensn::log.LogAllVerbose2() << "MeshHandlerCreate: Mesh Handler " << index << " created\n";

  return 1;
}

int
chiMeshHandlerSetCurrent(lua_State* L)
{
  int num_args = lua_gettop(L);
  if (num_args != 1) LuaPostArgAmountError("chiMeshHandlerSetCurrent", 1, num_args);

  int handle = lua_tonumber(L, 1);

  if ((handle < 0) or (handle >= opensn::meshhandler_stack.size()))
  {
    opensn::log.LogAllError() << "Invalid handle to mesh handler specified "
                              << "in call to chiMeshHandlerSetCurrent";
    opensn::Exit(EXIT_FAILURE);
  }

  opensn::current_mesh_handler = handle;

  opensn::log.LogAllVerbose2() << "chiMeshHandlerSetCurrent: set to " << handle;

  return 0;
}
