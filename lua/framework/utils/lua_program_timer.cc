// SPDX-FileCopyrightText: 2024 The OpenSn Authors <https://open-sn.github.io/opensn/>
// SPDX-License-Identifier: MIT

#include "framework/console/console.h"

#include "framework/lua.h"

#include "framework/runtime.h"
#include "framework/utils/timer.h"

using namespace opensn;

namespace opensnlua
{

/**Returns the program time as determined from the home location (involves a
 * collective broadcast).*/
int ProgramTime(lua_State* L);

RegisterLuaFunction(ProgramTime);

int
ProgramTime(lua_State* L)
{
  double time;
  if (opensn::mpi_comm.rank() == 0)
    time = opensn::program_timer.GetTime() / 1000.0;

  opensn::mpi_comm.broadcast(time, 0);

  return LuaReturn(L, time);
}

} // namespace opensnlua
