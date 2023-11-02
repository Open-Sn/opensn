#ifndef CHITECH_VOLUMEMESHER_LUA_H
#define CHITECH_VOLUMEMESHER_LUA_H

#include "chi_lua.h"

int chiVolumeMesherCreate(lua_State *L);
int chiVolumeMesherExecute(lua_State *L);
int chiVolumeMesherSetProperty(lua_State *L);

int chiVolumeMesherSetKBAPartitioningPxPyPz(lua_State *L);
int chiVolumeMesherSetKBACutsX(lua_State *L);
int chiVolumeMesherSetKBACutsY(lua_State *L);
int chiVolumeMesherSetKBACutsZ(lua_State *L);

int chiVolumeMesherSetMatIDToAll(lua_State* L);
int chiVolumeMesherSetupOrthogonalBoundaries(lua_State* L);

#endif //CHITECH_VOLUMEMESHER_LUA_H
