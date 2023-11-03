#pragma once

#include "chi_lua.h"

int chiGetFieldFunctionHandleByName(lua_State* L);
int chiExportFieldFunctionToVTK(lua_State* L);
int chiExportMultiFieldFunctionToVTK(lua_State* L);


