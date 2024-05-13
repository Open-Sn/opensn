// SPDX-FileCopyrightText: 2024 The OpenSn Authors <https://open-sn.github.io/opensn/>
// SPDX-License-Identifier: MIT

#include "framework/lua.h"
#include "framework/mesh/unpartitioned_mesh/unpartitioned_mesh.h"
#include "framework/runtime.h"
#include "framework/logging/log.h"
#include "framework/mesh/io/mesh_io.h"
#include "unpartition_mesh_lua_utils.h"
#include "framework/console/console.h"

using namespace opensn;

namespace opensnlua
{

RegisterLuaFunctionInNamespace(MeshUnpartitionedMeshFromVTU, mesh, UnpartitionedMeshFromVTU);
RegisterLuaFunctionInNamespace(MeshUnpartitionedMeshFromPVTU, mesh, UnpartitionedMeshFromPVTU);
RegisterLuaFunctionInNamespace(MeshUnpartitionedMeshFromEnsightGold,
                               mesh,
                               UnpartitionedMeshFromEnsightGold);
RegisterLuaFunctionInNamespace(MeshUnpartitionedMeshFromWavefrontOBJ,
                               mesh,
                               UnpartitionedMeshFromWavefrontOBJ);
RegisterLuaFunctionInNamespace(MeshUnpartitionedMeshFromMshFormat,
                               mesh,
                               UnpartitionedMeshFromMshFormat);
RegisterLuaFunctionInNamespace(MeshUnpartitionedMeshFromExodusII,
                               mesh,
                               UnpartitionedMeshFromExodusII);

int
MeshUnpartitionedMeshFromVTU(lua_State* L)
{
  const std::string fname = "mesh.UnpartitionedMeshFromVTU";
  LuaCheckArgs<std::string>(L, fname);

  auto file_name = LuaArg<std::string>(L, 1);
  auto field = LuaArgOptional<std::string>(L, 2, std::string(""));

  UnpartitionedMesh::Options options;
  options.file_name = file_name;
  options.material_id_fieldname = field;
  options.boundary_id_fieldname = field;

  auto new_object = MeshIO::FromVTU(options);
  opensn::unpartitionedmesh_stack.emplace_back(new_object);

  auto index = opensn::unpartitionedmesh_stack.size() - 1;
  return LuaReturn(L, index);
}

int
MeshUnpartitionedMeshFromPVTU(lua_State* L)
{
  const std::string func_name = "mesh.UnpartitionedMeshFromPVTU";
  LuaCheckArgs<std::string>(L, func_name);

  auto file_name = LuaArg<std::string>(L, 1);
  auto field = LuaArgOptional<std::string>(L, 2, std::string(""));

  opensn::UnpartitionedMesh::Options options;
  options.file_name = file_name;
  options.material_id_fieldname = field;
  options.boundary_id_fieldname = field;

  auto new_object = MeshIO::FromPVTU(options);
  opensn::unpartitionedmesh_stack.emplace_back(new_object);

  auto index = opensn::unpartitionedmesh_stack.size() - 1;
  return LuaReturn(L, index);
}

int
MeshUnpartitionedMeshFromEnsightGold(lua_State* L)
{
  const std::string func_name = "mesh.UnpartitionedMeshFromEnsightGold";
  LuaCheckArgs<std::string>(L, func_name);

  auto file_name = LuaArg<std::string>(L, 1);
  auto scale = LuaArgOptional<double>(L, 2, 1.0);

  opensn::UnpartitionedMesh::Options options;
  options.file_name = file_name;
  options.scale = scale;

  auto new_object = MeshIO::FromEnsightGold(options);
  opensn::unpartitionedmesh_stack.emplace_back(new_object);

  auto index = opensn::unpartitionedmesh_stack.size() - 1;
  return LuaReturn(L, index);
}

int
MeshUnpartitionedMeshFromWavefrontOBJ(lua_State* L)
{
  const std::string func_name = "mesh.UnpartitionedMeshFromWavefrontOBJ";
  LuaCheckArgs<std::string>(L, func_name);

  auto file_name = LuaArg<std::string>(L, 1);

  opensn::UnpartitionedMesh::Options options;
  options.file_name = file_name;

  auto new_object = MeshIO::FromOBJ(options);
  opensn::unpartitionedmesh_stack.emplace_back(new_object);

  auto index = opensn::unpartitionedmesh_stack.size() - 1;
  return LuaReturn(L, index);
}

int
MeshUnpartitionedMeshFromMshFormat(lua_State* L)
{
  const std::string func_name = "mesh.UnpartitionedMeshFromMshFormat";
  LuaCheckArgs<std::string>(L, func_name);

  auto file_name = LuaArg<std::string>(L, 1);

  opensn::UnpartitionedMesh::Options options;
  options.file_name = file_name;

  auto new_object = MeshIO::FromGmsh(options);
  opensn::unpartitionedmesh_stack.emplace_back(new_object);

  auto index = opensn::unpartitionedmesh_stack.size() - 1;
  return LuaReturn(L, index);
}

int
MeshUnpartitionedMeshFromExodusII(lua_State* L)
{
  const std::string func_name = "mesh.UnpartitionedMeshFromExodusII";
  LuaCheckArgs<std::string>(L, func_name);

  auto file_name = LuaArg<std::string>(L, 1);
  auto scale = LuaArgOptional<double>(L, 2, 1.0);

  opensn::UnpartitionedMesh::Options options;
  options.file_name = file_name;
  options.scale = scale;

  auto new_object = MeshIO::FromExodusII(options);
  opensn::unpartitionedmesh_stack.emplace_back(new_object);

  auto index = opensn::unpartitionedmesh_stack.size() - 1;
  return LuaReturn(L, index);
}

} // namespace opensnlua
