#include "framework/lua.h"

#include "framework/mesh/mesh_handler/mesh_handler.h"
#include "framework/mesh/mesh_continuum/mesh_continuum.h"

int
MeshHandlerExportMeshToObj(lua_State* L)
{
  // Check arguments
  const std::string fname = __FUNCTION__;
  const int num_args = lua_gettop(L);
  if (num_args < 1) LuaPostArgAmountError(fname, 1, num_args);

  const std::string file_name = lua_tostring(L, 1);

  bool per_material = false;
  if (num_args == 2) per_material = lua_toboolean(L, 2);

  // Get current handler
  auto& cur_hndlr = opensn::GetCurrentHandler();

  auto& grid = cur_hndlr.GetGrid();
  grid->ExportCellsToObj(file_name.c_str(), per_material);

  return 0;
}

int
MeshHandlerExportMeshToVTK(lua_State* L)
{
  // Check arguments
  const std::string fname = __FUNCTION__;
  const int num_args = lua_gettop(L);
  if (num_args != 1) LuaPostArgAmountError(fname, 1, num_args);

  const std::string file_name = lua_tostring(L, 1);

  // Get current handler
  auto& cur_hndlr = opensn::GetCurrentHandler();

  auto& grid = cur_hndlr.GetGrid();
  grid->ExportCellsToVTK(file_name);

  return 0;
}

int
MeshHandlerExportMeshToExodus(lua_State* L)
{
  // Check arguments
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

  // Get current handler
  auto& cur_hndlr = opensn::GetCurrentHandler();

  auto& grid = cur_hndlr.GetGrid();
  grid->ExportCellsToExodus(file_name, suppress_nodesets, suppress_sidesets);

  return 0;
}
