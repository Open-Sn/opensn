// SPDX-FileCopyrightText: 2024 The OpenSn Authors <https://open-sn.github.io/opensn/>
// SPDX-License-Identifier: MIT

#include "lua/framework/exceptions/invalid_lua_argument.h"
#include "lua/framework/lua.h"
#include <exception>

namespace opensnlua
{

InvalidLuaArgument::InvalidLuaArgument(lua_State* L, const std::string& what) : std::exception()
{
  auto dbg = LuaDebug(L);
  text_ = std::string(dbg.short_src) + ":" + std::to_string(dbg.currentline) + ": " + what;
}

InvalidLuaArgument::InvalidLuaArgument(lua_State* L, const char* what) : std::exception()
{
  auto dbg = LuaDebug(L);
  text_ =
    std::string(dbg.short_src) + ":" + std::to_string(dbg.currentline) + ": " + std::string(what);
}

const char*
InvalidLuaArgument::what() const noexcept
{
  return text_.c_str();
}

} // namespace opensnlua
