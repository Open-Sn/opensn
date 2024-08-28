// SPDX-FileCopyrightText: 2024 The OpenSn Authors <https://open-sn.github.io/opensn/>
// SPDX-License-Identifier: MIT

#include "lua/framework/lua.h"
#include "lua/framework/mpi/mpi_barrier.h"
#include "lua/framework/console/console.h"
#include "framework/runtime.h"

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
