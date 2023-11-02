#include "chi_lua.h"

#include "../chi_surfacemesh.h"
#include "mesh/MeshHandler/chi_meshhandler.h"

#include "chi_runtime.h"

#include "chi_log.h"
#include "lua_surface_mesh.h"
#include "console/chi_console.h"

RegisterLuaFunctionAsIs(chiSurfaceMeshCreate);

/** \defgroup LuaSurfaceMesh Surface Meshes
 * \ingroup LuaMesh
 */

//#############################################################################
/** Creates a new empty surface mesh.

### Example
Example usage:
\code
surfmesh = chiSurfaceMeshCreate()
\endcode

\return Handle int Handle to the created surface mesh.
\ingroup LuaSurfaceMesh
\author Jan*/
int
chiSurfaceMeshCreate(lua_State* L)
{
  auto new_mesh = new chi_mesh::SurfaceMesh;

  Chi::surface_mesh_stack.emplace_back(new_mesh);

  size_t index = Chi::surface_mesh_stack.size() - 1;
  lua_pushnumber(L, static_cast<lua_Number>(index));

  Chi::log.LogAllVerbose2() << "chiSurfaceMeshCreate: "
                               "Empty SurfaceMesh object, "
                            << index << ", created" << std::endl;

  return 1;
}