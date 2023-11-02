#ifndef CHITECH_FFINTERPOL_LUA_H
#define CHITECH_FFINTERPOL_LUA_H

#include "chi_lua.h"

int chiFFInterpolationCreate(lua_State* L);
int chiFFInterpolationSetProperty(lua_State* L);
int chiFFInterpolationInitialize(lua_State* L);
int chiFFInterpolationExecute(lua_State* L);
int chiFFInterpolationExportPython(lua_State* L);
int chiFFInterpolationGetValue(lua_State* L);

#endif // CHITECH_FFINTERPOL_LUA_H
