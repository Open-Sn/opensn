#ifndef CHI_PHYSICS_LUA_UTILS_H
#define CHI_PHYSICS_LUA_UTILS_H

#include "chi_lua.h"

int chiPhysicsAddMaterial(lua_State* L);
int chiPhysicsMaterialAddProperty(lua_State* L);
int chiPhysicsMaterialSetProperty(lua_State* L);
int chiPhysicsMaterialGetProperty(lua_State* L);

#endif // CHI_PHYSICS_LUA_UTILS_H