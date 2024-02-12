#include "framework/lua.h"

#include "framework/mesh/mesh_handler/mesh_handler.h"
#include "framework/runtime.h"

#include "framework/logging/log.h"

#include <iostream>
#include "mesh_handler_lua.h"
#include "framework/console/console.h"

using namespace opensn;

RegisterLuaFunctionAsIs(MeshHandlerExportMeshToObj);
RegisterLuaFunctionAsIs(MeshHandlerExportMeshToVTK);
RegisterLuaFunctionAsIs(MeshHandlerExportMeshToExodus);
