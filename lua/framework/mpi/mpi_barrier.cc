#include "opensn/framework/chi_lua.h"

#include "opensn/framework/chi_runtime.h"
#include "chi_mpi_lua.h"
#include "opensn/framework/console/chi_console.h"

namespace chi_mpi_utils
{

RegisterLuaFunctionAsIs(chiMPIBarrier);

// #############################################################################
/** Blocks until all processes in the communicator have reached this routine.

\ingroup chiMPI
\author Jan*/
int
chiMPIBarrier(lua_State* L)
{

  MPI_Barrier(Chi::mpi.comm);
  return 0;
}

} // namespace chi_mpi_utils
