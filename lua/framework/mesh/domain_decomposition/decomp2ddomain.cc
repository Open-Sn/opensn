#include "framework/chi_lua.h"
#include "framework/chi_runtime.h"

#include "framework/mesh/SurfaceMesh/chi_surfacemesh.h"
#include "domaindecomp_lua.h"
#include "framework/console/chi_console.h"

RegisterLuaFunctionAsIs(chiDecomposeSurfaceMeshPxPy);

/** \defgroup LuaDomainDecomposition Domain decomposition
 * \ingroup LuaMesh
 */

/** Decomposes a surface mesh into block px py elements.
 * \image html "InProgressImage.png" width=200px
 *
\param Surface mesh handler
\param Px int Number of divisions in x.
\param Py int Number of divisions in y.

\ingroup LuaDomainDecomposition
\author Jan*/
int
chiDecomposeSurfaceMeshPxPy(lua_State* L)
{
  int num_args = lua_gettop(L);

  if (num_args != 3) LuaPostArgAmountError("chiDecomposeSurfaceMeshPxPy", 3, num_args);

  //================================================== Extract arguments
  int surface_hndl = lua_tonumber(L, 1);
  int px = lua_tonumber(L, 2);
  int py = lua_tonumber(L, 3);

  auto& surf_mesh = Chi::GetStackItem<chi_mesh::SurfaceMesh>(Chi::surface_mesh_stack, surface_hndl);

  chi_mesh::DecomposeSurfaceMeshPxPy(surf_mesh, px, py);

  return 0;
}
