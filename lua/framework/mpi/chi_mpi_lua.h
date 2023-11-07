#pragma once

#include "framework/lua.h"

namespace chi_mpi_utils
{

/** Blocks until all processes in the communicator have reached this routine.
 *
 * \ingroup chiMPI
 * \author Jan
 */
int chiMPIBarrier(lua_State* L);

} // namespace chi_mpi_utils
