// SPDX-FileCopyrightText: 2024 The OpenSn Authors <https://open-sn.github.io/opensn/>
// SPDX-License-Identifier: MIT

#include "lua/framework/lua.h"
#include "framework/mesh/mesh.h"
#include "framework/mesh/mesh_continuum/mesh_continuum.h"
#include "framework/mesh/io/mesh_io.h"
#include "lua/framework/mesh/export/mesh_export.h"
#include "lua/framework/console/console.h"

namespace opensnlua
{

RegisterLuaFunctionInNamespace(MeshExportToOBJ, mesh, ExportToOBJ);
RegisterLuaFunctionInNamespace(MeshExportToPVTU, mesh, ExportToPVTU);
RegisterLuaFunctionInNamespace(MeshExportToExodusII, mesh, ExportToExodusII);

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
MeshExportToPVTU(lua_State* L)
{
  const std::string fname = "mesh.ExportToPVTU";
  LuaCheckArgs<std::string>(L, fname);

  const auto file_name = LuaArg<std::string>(L, 1);

  auto grid = opensn::GetCurrentMesh();
  opensn::MeshIO::ToPVTU(grid, file_name);

  return LuaReturn(L);
}

int
MeshExportToExodusII(lua_State* L)
{
  const std::string fname = "mesh.ExportToExodusII";
  LuaCheckArgs<std::string>(L, fname);

  const auto file_name = LuaArg<std::string>(L, 1);
  bool write_node_sets = LuaArgOptional<bool>(L, 2, true);
  bool write_side_sets = LuaArgOptional<bool>(L, 3, true);

  auto grid = opensn::GetCurrentMesh();
  opensn::MeshIO::ToExodusII(grid, file_name, write_node_sets, write_side_sets);

  return LuaReturn(L);
}

} // namespace opensnlua
