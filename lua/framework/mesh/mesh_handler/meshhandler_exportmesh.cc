#include "opensn/framework/chi_lua.h"

#include "opensn/framework/mesh/MeshHandler/chi_meshhandler.h"
#include "opensn/framework/mesh/MeshContinuum/chi_meshcontinuum.h"

//###################################################################
/**Exports the mesh to a wavefront.obj format.
\param FileName char Base name of the file to be used.
\param ExportByMaterial bool Default: False. Flag indicating whether to export
                     the extruder's surface mesh by material.
\ingroup LuaMeshHandler
*/
int
chiMeshHandlerExportMeshToObj(lua_State* L)
{
  //============================================= Check arguments
  const std::string fname = __FUNCTION__;
  const int num_args = lua_gettop(L);
  if (num_args < 1) LuaPostArgAmountError(fname, 1, num_args);

  const std::string file_name = lua_tostring(L, 1);

  bool per_material = false;
  if (num_args == 2) per_material = lua_toboolean(L, 2);

  //============================================= Get current handler
  auto& cur_hndlr = chi_mesh::GetCurrentHandler();

  auto& grid = cur_hndlr.GetGrid();
  grid->ExportCellsToObj(file_name.c_str(), per_material);

  return 0;
}

//###################################################################
/**Exports the mesh to vtu format.
\param FileName char Base name of the file to be used.
\ingroup LuaMeshHandler
*/
int
chiMeshHandlerExportMeshToVTK(lua_State* L)
{
  //============================================= Check arguments
  const std::string fname = __FUNCTION__;
  const int num_args = lua_gettop(L);
  if (num_args != 1) LuaPostArgAmountError(fname, 1, num_args);

  const std::string file_name = lua_tostring(L, 1);

  //============================================= Get current handler
  auto& cur_hndlr = chi_mesh::GetCurrentHandler();

  auto& grid = cur_hndlr.GetGrid();
  grid->ExportCellsToVTK(file_name);

  return 0;
}

//###################################################################
/**Exports the mesh to exodus format (.e extensions).
\param FileName char Base name of the file to be used.
\param suppress_nodesets bool Optional. Flag to suppress exporting nodesets.
                              Default = `false`.
\param suppress_sidesets bool Optional. Flag to suppress exporting sidesets.
                              Default = `false`.
\ingroup LuaMeshHandler
*/
int
chiMeshHandlerExportMeshToExodus(lua_State* L)
{
  //============================================= Check arguments
  const std::string fname = __FUNCTION__;
  const int num_args = lua_gettop(L);
  if (num_args < 1) LuaPostArgAmountError(fname, 1, num_args);

  const std::string file_name = lua_tostring(L, 1);

  bool suppress_nodesets = false;
  bool suppress_sidesets = false;
  if (num_args >= 2)
  {
    LuaCheckBoolValue(fname, L, 2);
    suppress_nodesets = lua_toboolean(L, 2);
  }

  if (num_args == 3)
  {
    LuaCheckBoolValue(fname, L, 3);
    suppress_sidesets = lua_toboolean(L, 3);
  }

  //============================================= Get current handler
  auto& cur_hndlr = chi_mesh::GetCurrentHandler();

  auto& grid = cur_hndlr.GetGrid();
  grid->ExportCellsToExodus(file_name, suppress_nodesets, suppress_sidesets);

  return 0;
}
