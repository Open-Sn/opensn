#include "framework/lua.h"

#include "framework/runtime.h"
#include "mpi_lua.h"
#include "framework/console/console.h"

namespace opensnlua
{

RegisterLuaFunctionAsIs(chiMPIBarrier);

int
chiMPIBarrier(lua_State* L)
{

  opensn::mpi_comm.barrier();
  return 0;
}

} // namespace opensnlua
