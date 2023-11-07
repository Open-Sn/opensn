#include "framework/lua.h"
#include "framework/runtime.h"

#include "framework/mesh/surface_mesh/surface_mesh.h"
#include "domaindecomp_lua.h"
#include "framework/console/console.h"

RegisterLuaFunctionAsIs(chiDecomposeSurfaceMeshPxPy);

int
chiDecomposeSurfaceMeshPxPy(lua_State* L)
{
  int num_args = lua_gettop(L);

  if (num_args != 3) LuaPostArgAmountError("chiDecomposeSurfaceMeshPxPy", 3, num_args);

  // Extract arguments
  int surface_hndl = lua_tonumber(L, 1);
  int px = lua_tonumber(L, 2);
  int py = lua_tonumber(L, 3);

  auto& surf_mesh = Chi::GetStackItem<chi_mesh::SurfaceMesh>(Chi::surface_mesh_stack, surface_hndl);

  chi_mesh::DecomposeSurfaceMeshPxPy(surf_mesh, px, py);

  return 0;
}
