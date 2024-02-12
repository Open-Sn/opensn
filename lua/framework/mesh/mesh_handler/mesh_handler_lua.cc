#include "framework/lua.h"

#include "framework/mesh/mesh_handler/mesh_handler.h"
#include "framework/runtime.h"

#include "framework/logging/log.h"

#include <iostream>
#include "mesh_handler_lua.h"
#include "framework/console/console.h"

using namespace opensn;

RegisterLuaFunctionAsIs(MeshHandlerCreate);
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
