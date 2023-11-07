#include "framework/lua.h"

#include "framework/runtime.h"
#include "chi_mpi_lua.h"
#include "framework/console/console.h"

namespace chi_mpi_utils
{

RegisterLuaFunctionAsIs(chiMPIBarrier);

int
chiMPIBarrier(lua_State* L)
{

  MPI_Barrier(Chi::mpi.comm);
  return 0;
}

} // namespace chi_mpi_utils
