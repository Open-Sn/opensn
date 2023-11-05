#pragma once

#include "chi_lua.h"

int chiPhysicsTransportXSCreate(lua_State* L);
int chiPhysicsTransportXSSet(lua_State* L);
int chiPhysicsTransportXSMakeCombined(lua_State* L);
int chiPhysicsTransportXSSetCombined(lua_State* L);
int chiPhysicsTransportXSGet(lua_State* L);
int chiPhysicsTransportXSExportToChiTechFormat(lua_State* L);
