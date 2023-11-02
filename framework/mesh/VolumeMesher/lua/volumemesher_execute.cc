#include "chi_lua.h"

#include "mesh/MeshHandler/chi_meshhandler.h"
#include "mesh/VolumeMesher/chi_volumemesher.h"

#include "chi_runtime.h"
#include "chi_log.h"

#include "utils/chi_timer.h"

#include "volumemesher_lua.h"
#include "console/chi_console.h"

#include <iomanip>
#include <iostream>

RegisterLuaFunctionAsIs(chiVolumeMesherExecute);

//#############################################################################
/** Executes the volume meshing pipeline.

\ingroup LuaVolumeMesher
\author Jan*/
int chiVolumeMesherExecute(lua_State *L)
{
  auto& cur_hndlr = chi_mesh::GetCurrentHandler();

  //Get memory before
  chi::CSTMemory mem_before = chi::Console::GetMemoryUsage();

  cur_hndlr.GetVolumeMesher().Execute();

  //Get memory usage
  chi::CSTMemory mem_after = chi::Console::GetMemoryUsage();

  std::stringstream mem_string;
  mem_string
  << " Memory used = " << std::setprecision(3)
  << mem_after.memory_mbytes - mem_before.memory_mbytes
  << " MB\n"
  << "Total process memory used after meshing "
  << mem_after.memory_mbytes
  << " MB";

  Chi::log.Log()
    << Chi::program_timer.GetTimeString()
    << " chiVolumeMesherExecute: Volume meshing completed."
    << mem_string.str()
    << std::endl;

  return 0;
}