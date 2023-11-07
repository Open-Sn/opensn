#include "framework/lua.h"

#include "framework/mesh/mesh_handler/mesh_handler.h"
#include "framework/mesh/volume_mesher/volume_mesher.h"

#include "framework/runtime.h"
#include "framework/logging/log.h"

#include "framework/utils/timer.h"

#include "volume_mesher_lua.h"
#include "framework/console/console.h"

#include <iomanip>
#include <iostream>

RegisterLuaFunctionAsIs(chiVolumeMesherExecute);

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
