#include "framework/chi_lua.h"

#include "framework/mesh/SurfaceMesher/surfacemesher.h"
#include "framework/mesh/MeshHandler/chi_meshhandler.h"

#include "framework/chi_runtime.h"
#include "framework/logging/chi_log.h"
#include "surfmesher_lua.h"
#include "framework/console/chi_console.h"

RegisterLuaFunctionAsIs(chiSurfaceMesherExecute);

int
chiSurfaceMesherExecute(lua_State* L)
{
  auto& cur_hndlr = chi_mesh::GetCurrentHandler();
  Chi::log.LogAllVerbose2() << "Executing surface mesher\n";

  cur_hndlr.GetSurfaceMesher().Execute();

  Chi::log.LogAllVerbose2() << "chiSurfaceMesherExecute: Surface mesher execution completed."
                            << std::endl;

  return 0;
}
