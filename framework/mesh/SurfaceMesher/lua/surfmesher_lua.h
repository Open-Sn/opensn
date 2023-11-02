#ifndef CHITECH_SURFMESHER_LUA_H
#define CHITECH_SURFMESHER_LUA_H

#include "chi_lua.h"

int chiSurfaceMesherCreate(lua_State* L);
int chiSurfaceMesherExecute(lua_State* L);
int chiSurfaceMesherSetProperty(lua_State* L);
// int chiSurfaceMesherExportToObj(lua_State* L);

#endif // CHITECH_SURFMESHER_LUA_H
