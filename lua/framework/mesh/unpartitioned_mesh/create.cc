#include "framework/lua.h"
#include "framework/memory_usage.h"
#include "framework/mesh/unpartitioned_mesh/unpartitioned_mesh.h"
#include "framework/runtime.h"
#include "framework/logging/log.h"

#include "unpartition_mesh_lua_utils.h"
#include "framework/console/console.h"

using namespace opensn;

namespace opensnlua
{

RegisterLuaFunctionNamespace(MeshUnpartitionedMeshFromVTU, mesh, UnpartitionedMeshFromVTU);
RegisterLuaFunctionNamespace(MeshUnpartitionedMeshFromPVTU, mesh, UnpartitionedMeshFromPVTU);
RegisterLuaFunctionNamespace(MeshUnpartitionedMeshFromEnsightGold,
                             mesh,
                             UnpartitionedMeshFromEnsightGold);
RegisterLuaFunctionNamespace(MeshUnpartitionedMeshFromWavefrontOBJ,
                             mesh,
                             UnpartitionedMeshFromWavefrontOBJ);
RegisterLuaFunctionNamespace(MeshUnpartitionedMeshFromMshFormat,
                             mesh,
                             UnpartitionedMeshFromMshFormat);
RegisterLuaFunctionNamespace(MeshUnpartitionedMeshFromExodusII,
                             mesh,
                             UnpartitionedMeshFromExodusII);

int
MeshUnpartitionedMeshFromVTU(lua_State* L)
{
  const std::string func_name = __FUNCTION__;
  int num_args = lua_gettop(L);
  if (num_args < 1)
    LuaPostArgAmountError(func_name, 1, num_args);

  auto file_name = LuaArg<std::string>(L, 1);
  auto field = LuaArgOptional<std::string>(L, 2, std::string(""));
  auto new_object = new UnpartitionedMesh;

  UnpartitionedMesh::Options options;
  options.file_name = file_name;
  options.material_id_fieldname = field;
  options.boundary_id_fieldname = field;

  new_object->ReadFromVTU(options);

  opensn::unpartitionedmesh_stack.emplace_back(new_object);

  lua_pushinteger(L, static_cast<lua_Integer>(opensn::unpartitionedmesh_stack.size() - 1));

  return 1;
}

int
MeshUnpartitionedMeshFromPVTU(lua_State* L)
{
  const std::string func_name = __FUNCTION__;
  int num_args = lua_gettop(L);
  if (num_args < 1)
    LuaPostArgAmountError(func_name, 1, num_args);

  auto file_name = LuaArg<std::string>(L, 1);
  auto field = LuaArgOptional<std::string>(L, 2, std::string(""));
  auto new_object = new opensn::UnpartitionedMesh;

  opensn::UnpartitionedMesh::Options options;
  options.file_name = file_name;
  options.material_id_fieldname = field;
  options.boundary_id_fieldname = field;

  new_object->ReadFromPVTU(options);

  opensn::unpartitionedmesh_stack.emplace_back(new_object);

  lua_pushinteger(L, static_cast<lua_Integer>(opensn::unpartitionedmesh_stack.size() - 1));

  return 1;
}

int
MeshUnpartitionedMeshFromEnsightGold(lua_State* L)
{
  const std::string func_name = __FUNCTION__;
  int num_args = lua_gettop(L);
  if (num_args < 1)
    LuaPostArgAmountError(func_name, 1, num_args);

  auto file_name = LuaArg<std::string>(L, 1);
  auto scale = LuaArgOptional<double>(L, 2, 1.0);
  auto new_object = new opensn::UnpartitionedMesh;

  opensn::UnpartitionedMesh::Options options;
  options.file_name = file_name;
  options.scale = scale;

  new_object->ReadFromEnsightGold(options);

  opensn::unpartitionedmesh_stack.emplace_back(new_object);

  lua_pushinteger(L, static_cast<lua_Integer>(opensn::unpartitionedmesh_stack.size() - 1));

  return 1;
}

int
MeshUnpartitionedMeshFromWavefrontOBJ(lua_State* L)
{
  const std::string func_name = __FUNCTION__;
  int num_args = lua_gettop(L);
  if (num_args < 1)
    LuaPostArgAmountError(func_name, 1, num_args);

  auto file_name = LuaArg<std::string>(L, 1);

  auto new_object = new opensn::UnpartitionedMesh;

  opensn::UnpartitionedMesh::Options options;
  options.file_name = file_name;

  new_object->ReadFromWavefrontOBJ(options);

  opensn::unpartitionedmesh_stack.emplace_back(new_object);

  lua_pushinteger(L, static_cast<lua_Integer>(opensn::unpartitionedmesh_stack.size() - 1));

  return 1;
}

int
MeshUnpartitionedMeshFromMshFormat(lua_State* L)
{
  const std::string func_name = __FUNCTION__;
  int num_args = lua_gettop(L);
  if (num_args < 1)
    LuaPostArgAmountError(func_name, 1, num_args);

  auto file_name = LuaArg<std::string>(L, 1);

  auto new_object = new opensn::UnpartitionedMesh;

  opensn::UnpartitionedMesh::Options options;
  options.file_name = file_name;

  new_object->ReadFromMsh(options);

  opensn::unpartitionedmesh_stack.emplace_back(new_object);

  lua_pushinteger(L, static_cast<lua_Integer>(opensn::unpartitionedmesh_stack.size() - 1));

  return 1;
}

int
MeshUnpartitionedMeshFromExodusII(lua_State* L)
{
  const std::string func_name = __FUNCTION__;
  int num_args = lua_gettop(L);
  if (num_args < 1)
    LuaPostArgAmountError(func_name, 1, num_args);

  auto file_name = LuaArg<std::string>(L, 1);
  auto scale = LuaArgOptional<double>(L, 2, 1.0);
  auto new_object = new opensn::UnpartitionedMesh;

  opensn::UnpartitionedMesh::Options options;
  options.file_name = file_name;
  options.scale = scale;

  new_object->ReadFromExodus(options);

  opensn::unpartitionedmesh_stack.emplace_back(new_object);

  lua_pushinteger(L, static_cast<lua_Integer>(opensn::unpartitionedmesh_stack.size() - 1));

  return 1;
}

} // namespace opensnlua
