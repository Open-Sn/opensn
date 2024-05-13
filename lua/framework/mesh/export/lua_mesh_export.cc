// SPDX-FileCopyrightText: 2024 The OpenSn Authors <https://open-sn.github.io/opensn/>
// SPDX-License-Identifier: MIT

#include "framework/lua.h"
#include "framework/mesh/mesh.h"
#include "framework/mesh/mesh_continuum/mesh_continuum.h"
#include "framework/mesh/io/mesh_io.h"
#include "lua/framework/mesh/export/lua_mesh_export.h"
#include "lua/framework/console/console.h"

namespace opensnlua
{

RegisterLuaFunctionInNamespace(MeshExportToOBJ, mesh, ExportToOBJ);
RegisterLuaFunctionInNamespace(MeshExportToVTK, mesh, ExportToVTK);
RegisterLuaFunctionInNamespace(MeshExportToExodus, mesh, ExportToExodus);

int
MeshExportToOBJ(lua_State* L)
{
  const std::string fname = "mesh.ExportToOBJ";
  LuaCheckArgs<int>(L, fname);

  const auto file_name = LuaArg<std::string>(L, 1);
  bool per_material = LuaArgOptional<bool>(L, 2, false);

  auto grid = opensn::GetCurrentMesh();
  opensn::MeshIO::ToOBJ(grid, file_name.c_str(), per_material);

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
