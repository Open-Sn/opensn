// SPDX-FileCopyrightText: 2024 The OpenSn Authors <https://open-sn.github.io/opensn/>
// SPDX-License-Identifier: MIT

#include "framework/lua.h"

#include "framework/runtime.h"
#include "mpi_lua.h"
#include "framework/console/console.h"

namespace opensnlua
{

RegisterLuaFunction(MPIBarrier);

int
MPIBarrier(lua_State* L)
{
  opensn::mpi_comm.barrier();
  return LuaReturn(L);
}

} // namespace opensnlua
