#pragma once

#include "opensn/framework/chi_lua.h"

int chiPhysicsAddMaterial(lua_State* L);
int chiPhysicsMaterialAddProperty(lua_State* L);
int chiPhysicsMaterialSetProperty(lua_State* L);
int chiPhysicsMaterialGetProperty(lua_State* L);
