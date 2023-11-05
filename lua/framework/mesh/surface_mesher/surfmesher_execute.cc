#include "opensn/framework/chi_lua.h"

#include "opensn/framework/mesh/SurfaceMesher/surfacemesher.h"
#include "opensn/framework/mesh/MeshHandler/chi_meshhandler.h"

#include "opensn/framework/chi_runtime.h"
#include "opensn/framework/logging/chi_log.h"
#include "surfmesher_lua.h"
#include "opensn/framework/console/chi_console.h"

RegisterLuaFunctionAsIs(chiSurfaceMesherExecute);

//#############################################################################
/** Executes the surface meshing pipeline.

\ingroup LuaSurfaceMesher
\author Jan*/
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
