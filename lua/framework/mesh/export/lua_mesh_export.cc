#include "framework/lua.h"

#include "framework/mesh/mesh.h"
#include "framework/mesh/mesh_continuum/mesh_continuum.h"
#include "lua/framework/mesh/export/lua_mesh_export.h"
#include "lua/framework/console/console.h"

RegisterLuaFunctionNamespace(MeshExportToObj, mesh, ExportToObj);
RegisterLuaFunctionNamespace(MeshExportToVTK, mesh, ExportToVTK);
RegisterLuaFunctionNamespace(MeshExportToExodus, mesh, ExportToExodus);

int
MeshExportToObj(lua_State* L)
{
  // Check arguments
  const std::string fname = __FUNCTION__;
  const int num_args = lua_gettop(L);
  if (num_args < 1)
    LuaPostArgAmountError(fname, 1, num_args);

  const std::string file_name = lua_tostring(L, 1);

  bool per_material = false;
  if (num_args == 2)
    per_material = lua_toboolean(L, 2);

  auto grid = opensn::GetCurrentMesh();
  grid->ExportCellsToObj(file_name.c_str(), per_material);

  return 0;
}

int
MeshExportToVTK(lua_State* L)
{
  // Check arguments
  const std::string fname = __FUNCTION__;
  const int num_args = lua_gettop(L);
  if (num_args != 1)
    LuaPostArgAmountError(fname, 1, num_args);

  const std::string file_name = lua_tostring(L, 1);

  auto grid = opensn::GetCurrentMesh();
  grid->ExportCellsToVTK(file_name);

  return 0;
}

int
MeshExportToExodus(lua_State* L)
{
  // Check arguments
  const std::string fname = __FUNCTION__;
  const int num_args = lua_gettop(L);
  if (num_args < 1)
    LuaPostArgAmountError(fname, 1, num_args);

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

  auto grid = opensn::GetCurrentMesh();
  grid->ExportCellsToExodus(file_name, suppress_nodesets, suppress_sidesets);

  return 0;
}
