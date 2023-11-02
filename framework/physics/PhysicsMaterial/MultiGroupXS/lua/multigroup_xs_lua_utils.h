#ifndef CHITECH_XSECTIONS_LUA_UTILS_H
#define CHITECH_XSECTIONS_LUA_UTILS_H

#include "chi_lua.h"

int chiPhysicsTransportXSCreate(lua_State* L);
int chiPhysicsTransportXSSet(lua_State* L);
int chiPhysicsTransportXSMakeCombined(lua_State* L);
int chiPhysicsTransportXSSetCombined(lua_State* L);
int chiPhysicsTransportXSGet(lua_State* L);
int chiPhysicsTransportXSExportToChiTechFormat(lua_State* L);

#endif // CHITECH_XSECTIONS_LUA_UTILS_H
