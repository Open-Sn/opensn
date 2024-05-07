// SPDX-FileCopyrightText: 2024 The OpenSn Authors <https://open-sn.github.io/opensn/>
// SPDX-License-Identifier: MIT

#include "framework/console/console.h"

#include "framework/lua.h"

#include "framework/utils/timer.h"

namespace opensnlua
{

/**Makes the program sleep for the specified time in milliseconds.
 * \param time int Time in milliseconds to sleep for.
 * */
int Sleep(lua_State* L);

RegisterLuaFunction(Sleep);

int
Sleep(lua_State* L)
{
  const std::string fname = "Sleep";
  LuaCheckArgs<int64_t>(L, fname);

  const auto time_to_sleep = LuaArg<int64_t>(L, 1);

  opensn::Sleep(std::chrono::milliseconds(time_to_sleep));

  return LuaReturn(L);
}

} // namespace opensnlua
