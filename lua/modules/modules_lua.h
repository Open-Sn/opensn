#pragma once

#include "lua/framework/lua.h"

namespace opensnlua
{
void LoadRegisteredLuaItems();
void RegisterLuaEntities(lua_State* L);
} // namespace opensnlua
