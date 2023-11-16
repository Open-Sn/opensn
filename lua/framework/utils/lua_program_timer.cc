#include "framework/console/console.h"

#include "framework/lua.h"

#include "framework/runtime.h"
#include "framework/utils/timer.h"

using namespace opensn;

namespace opensnlua
{

/**Returns the program time as determined from the home location (involves a
 * collective broadcast).*/
int chiProgramTime(lua_State* L);

RegisterLuaFunctionAsIs(chiProgramTime);

int
chiProgramTime(lua_State* L)
{
  double time;
  if (opensn::Chi::mpi.location_id == 0) time = opensn::Chi::program_timer.GetTime() / 1000.0;

  MPI_Bcast(&time, 1, MPI_DOUBLE, 0, opensn::Chi::mpi.comm);

  lua_pushnumber(L, time);
  return 1;
}

} // namespace opensnlua
