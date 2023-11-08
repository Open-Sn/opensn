#include "framework/chi_lua.h"

#include "framework/mesh/MeshHandler/chi_meshhandler.h"
#include "framework/mesh/VolumeMesher/chi_volumemesher.h"

#include "framework/chi_runtime.h"
#include "framework/logging/chi_log.h"

#include "framework/utils/chi_timer.h"

#include "volumemesher_lua.h"
#include "framework/console/chi_console.h"

#include <iomanip>
#include <iostream>

RegisterLuaFunctionAsIs(chiVolumeMesherExecute);

// #############################################################################
/** Executes the volume meshing pipeline.

\ingroup LuaVolumeMesher
\author Jan*/
int
chiVolumeMesherExecute(lua_State* L)
{
  auto& cur_hndlr = chi_mesh::GetCurrentHandler();

  // Get memory before
  chi::CSTMemory mem_before = Chi::GetMemoryUsage();

  cur_hndlr.GetVolumeMesher().Execute();

  // Get memory usage
  chi::CSTMemory mem_after = Chi::GetMemoryUsage();

  std::stringstream mem_string;
  mem_string << " Memory used = " << std::setprecision(3)
             << mem_after.memory_mbytes - mem_before.memory_mbytes << " MB\n"
             << "Total process memory used after meshing " << mem_after.memory_mbytes << " MB";

  Chi::log.Log() << Chi::program_timer.GetTimeString()
                 << " chiVolumeMesherExecute: Volume meshing completed." << mem_string.str()
                 << std::endl;

  return 0;
}
