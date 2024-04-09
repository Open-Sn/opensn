#include "framework/lua.h"
#include "framework/mesh/mesh.h"
#include "framework/mesh/mesh_continuum/mesh_continuum.h"
#include "lua/framework/mesh/export/lua_mesh_export.h"
#include "lua/framework/console/console.h"

namespace opensnlua
{

RegisterLuaFunctionNamespace(MeshExportToObj, mesh, ExportToObj);
RegisterLuaFunctionNamespace(MeshExportToVTK, mesh, ExportToVTK);
RegisterLuaFunctionNamespace(MeshExportToExodus, mesh, ExportToExodus);

int
MeshExportToObj(lua_State* L)
{
  const std::string fname = "mesh.ExportToObj";
  LuaCheckArgs<int>(L, fname);

  const auto file_name = LuaArg<std::string>(L, 1);
  bool per_material = LuaArgOptional<bool>(L, 2, false);

  auto grid = opensn::GetCurrentMesh();
  grid->ExportCellsToObj(file_name.c_str(), per_material);

  return LuaReturn(L);
}

int
MeshExportToVTK(lua_State* L)
{
  const std::string fname = "mesh.ExportToVTK";
  LuaCheckArgs<std::string>(L, fname);

  const auto file_name = LuaArg<std::string>(L, 1);

  auto grid = opensn::GetCurrentMesh();
  grid->ExportCellsToVTK(file_name);

  return LuaReturn(L);
}

int
MeshExportToExodus(lua_State* L)
{
  const std::string fname = "mesh.ExportToExodus";
  LuaCheckArgs<std::string>(L, fname);

  const auto file_name = LuaArg<std::string>(L, 1);
  bool suppress_nodesets = LuaArgOptional<bool>(L, 2, false);
  bool suppress_sidesets = LuaArgOptional<bool>(L, 3, false);

  auto grid = opensn::GetCurrentMesh();
  grid->ExportCellsToExodus(file_name, suppress_nodesets, suppress_sidesets);

  return LuaReturn(L);
}

} // namespace opensnlua
