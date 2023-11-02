#include "chi_lua.h"

#include <iostream>
#include "../chi_surfacemesh.h"
#include "mesh/MeshHandler/chi_meshhandler.h"

#include "chi_runtime.h"

#include "chi_log.h"
#include "lua_surface_mesh.h"
#include "console/chi_console.h"

RegisterLuaFunctionAsIs(chiSurfaceMeshExportToObj);
RegisterLuaFunctionAsIs(chiSurfaceMeshExportPolyFile);

//############################################################################# Create
/** Exports mesh as a .obj format.
 *
\param SurfaceHandle int Handle to the surface on which the operation is to be performed.
\param FileName char* Path to the file to be exported.

\ingroup LuaSurfaceMesh
\author Jan*/
int chiSurfaceMeshExportToObj(lua_State* L)
{
  auto& cur_hndlr = chi_mesh::GetCurrentHandler();

  //============================================= Get arguments
  int num_args = lua_gettop(L);
  if (num_args != 2)
    LuaPostArgAmountError("chiSurfaceMeshExportObj",2, num_args);

  int handle = lua_tonumber(L,1);

  size_t length = 0;
  const char* temp = lua_tolstring(L, 2, &length);

  auto& surface_mesh = Chi::GetStackItem<chi_mesh::SurfaceMesh>(
    Chi::surface_mesh_stack, handle, __FUNCTION__);

  surface_mesh.ExportToOBJFile(temp);
  
  return 0;
}

//############################################################################# Create
/** Exports mesh as a .poly format.
 *
\param SurfaceHandle int Handle to the surface on which the operation is to be performed.
\param FileName char* Path and basename to the file to be exported.

\ingroup LuaSurfaceMesh
\author Jan*/
int chiSurfaceMeshExportPolyFile(lua_State* L)
{
  auto& cur_hndlr = chi_mesh::GetCurrentHandler();

  //============================================= Get arguments
  int num_args = lua_gettop(L);
  if (num_args != 2)
    LuaPostArgAmountError("chiSurfaceMeshExportPolyFile",2, num_args);

  int handle = lua_tonumber(L,1);

  size_t length = 0;
  const char* temp = lua_tolstring(L, 2, &length);

  auto& surface_mesh = Chi::GetStackItem<chi_mesh::SurfaceMesh>(
    Chi::surface_mesh_stack, handle, __FUNCTION__);

  surface_mesh.ExportToPolyFile(temp);
  return 0;
}
