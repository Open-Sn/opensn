// SPDX-FileCopyrightText: 2024 The OpenSn Authors <https://open-sn.github.io/opensn/>
// SPDX-License-Identifier: MIT

#include "framework/console/console.h"

#include "framework/lua.h"

#include "framework/runtime.h"

namespace opensnlua
{

/**
 * Gracefully exits OpenSn.
 * \param return_code int Return code, defaults to 0 (Success).
 */
int Exit(lua_State* L);

RegisterLuaFunction(Exit);

int
Exit(lua_State* L)
{
  auto return_code = LuaArgOptional<int>(L, 1, EXIT_SUCCESS);
  opensn::Exit(return_code);
  return LuaReturn(L);
}

} // namespace opensnlua
