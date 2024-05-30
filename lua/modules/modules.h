// SPDX-FileCopyrightText: 2024 The OpenSn Authors <https://open-sn.github.io/opensn/>
// SPDX-License-Identifier: MIT

#pragma once

#include "lua/framework/lua.h"

namespace opensnlua
{
void LoadRegisteredLuaItems();
void RegisterLuaEntities(lua_State* L);
} // namespace opensnlua
