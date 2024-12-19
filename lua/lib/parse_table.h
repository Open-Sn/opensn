// SPDX-FileCopyrightText: 2024 The OpenSn Authors <https://open-sn.github.io/opensn/>
// SPDX-License-Identifier: MIT

#pragma once

extern "C"
{
#include <lua.h>
#include <lualib.h>
#include <lauxlib.h>
}

namespace opensn
{
class ParameterBlock;
}

namespace luabridge
{
struct Result;
}

namespace opensnlua
{

void ParseTableKeys(lua_State* L, int idx, opensn::ParameterBlock& block);

luabridge::Result PushParameterBlock(lua_State* L, const opensn::ParameterBlock& block, int level);

} // namespace opensnlua
