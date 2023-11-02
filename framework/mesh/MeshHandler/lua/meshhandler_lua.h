#ifndef CHITECH_MESHHANDLER_LUA_H
#define CHITECH_MESHHANDLER_LUA_H

#include "chi_lua.h"

int chiMeshHandlerCreate(lua_State* L);
int chiMeshHandlerSetCurrent(lua_State* L);
int chiMeshHandlerExportMeshToObj(lua_State* L);
int chiMeshHandlerExportMeshToVTK(lua_State* L);
int chiMeshHandlerExportMeshToExodus(lua_State* L);

#endif // CHITECH_MESHHANDLER_LUA_H
