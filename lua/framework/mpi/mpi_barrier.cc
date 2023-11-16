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

  MPI_Barrier(opensn::Chi::mpi.comm);
  return 0;
}

} // namespace opensnlua
