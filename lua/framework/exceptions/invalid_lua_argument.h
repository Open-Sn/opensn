// SPDX-FileCopyrightText: 2024 The OpenSn Authors <https://open-sn.github.io/opensn/>
// SPDX-License-Identifier: MIT

#pragma once

extern "C"
{
#include <lua.h>
}
#include <exception>
#include <string>

namespace opensnlua
{

class InvalidLuaArgument : public std::exception
{
public:
  InvalidLuaArgument(lua_State* L, const std::string& what);
  InvalidLuaArgument(lua_State* L, const char* what);

  const char* what() const noexcept override;

private:
  std::string text_;
};

} // namespace opensnlua
