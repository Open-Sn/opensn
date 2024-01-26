#include "framework/lua.h"

#include "framework/mesh/surface_mesh/surface_mesh.h"

#include "framework/runtime.h"

#include "framework/logging/log.h"

using namespace opensn;

/** Exports all open edges of a surface mesh to file. This is used mostly
 * for graphical error checking.
 *
 * \param SurfaceHandle int Handle to the surface on which the operation is to be
 * performed. \param FileName char Filename to which the edges are to be exported.
 *
 * \ingroup LuaSurfaceMesh
 * \author Jan
 */
int
chiSurfaceMeshExtractOpenEdgesToObj(lua_State* L)
{
  int num_args = lua_gettop(L);
  if (num_args != 2) LuaPostArgAmountError("chiSurfaceMeshExtractOpenEdgesToObj", 2, num_args);

  auto& cur_hndlr = opensn::GetCurrentHandler();

  int surf_handle = lua_tonumber(L, 1);
  const char* file_name = lua_tostring(L, 2);

  auto& surface_mesh =
    opensn::GetStackItem<SurfaceMesh>(opensn::surface_mesh_stack, surf_handle, __FUNCTION__);

  surface_mesh.ExtractOpenEdgesToObj(file_name);
  return 0;
}
