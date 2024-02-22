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

RegisterLuaFunctionAsIs(UnpartitionedMeshFromVTU);
RegisterLuaFunctionAsIs(UnpartitionedMeshFromPVTU);
RegisterLuaFunctionAsIs(UnpartitionedMeshFromEnsightGold);
RegisterLuaFunctionAsIs(UnpartitionedMeshFromWavefrontOBJ);
RegisterLuaFunctionAsIs(UnpartitionedMeshFromMshFormat);
RegisterLuaFunctionAsIs(UnpartitionedMeshFromExodusII);

int
UnpartitionedMeshFromVTU(lua_State* L)
{
  const std::string func_name = __FUNCTION__;
  int num_args = lua_gettop(L);
  if (num_args < 1)
    LuaPostArgAmountError(func_name, 1, num_args);

  LuaCheckNilValue(func_name, L, 1);
  if (num_args >= 2)
    LuaCheckNilValue(func_name, L, 2);

  const char* temp = lua_tostring(L, 1);
  const char* field = "";
  if (num_args >= 2)
    field = lua_tostring(L, 2);
  auto new_object = new UnpartitionedMesh;

  UnpartitionedMesh::Options options;
  options.file_name = std::string(temp);
  options.material_id_fieldname = field;
  options.boundary_id_fieldname = field;

  new_object->ReadFromVTU(options);

  opensn::unpartitionedmesh_stack.emplace_back(new_object);

  lua_pushnumber(L, static_cast<lua_Number>(opensn::unpartitionedmesh_stack.size() - 1));

  return 1;
}

int
UnpartitionedMeshFromPVTU(lua_State* L)
{
  const std::string func_name = __FUNCTION__;
  int num_args = lua_gettop(L);
  if (num_args < 1)
    LuaPostArgAmountError(func_name, 1, num_args);

  LuaCheckNilValue(func_name, L, 1);
  if (num_args >= 2)
    LuaCheckNilValue(func_name, L, 2);

  const char* temp = lua_tostring(L, 1);
  const char* field = "";
  if (num_args >= 2)
    field = lua_tostring(L, 2);
  auto new_object = new opensn::UnpartitionedMesh;

  opensn::UnpartitionedMesh::Options options;
  options.file_name = std::string(temp);
  options.material_id_fieldname = field;
  options.boundary_id_fieldname = field;

  new_object->ReadFromPVTU(options);

  opensn::unpartitionedmesh_stack.emplace_back(new_object);

  lua_pushnumber(L, static_cast<lua_Number>(opensn::unpartitionedmesh_stack.size() - 1));

  return 1;
}

int
UnpartitionedMeshFromEnsightGold(lua_State* L)
{
  const std::string func_name = __FUNCTION__;
  int num_args = lua_gettop(L);
  if (num_args < 1)
    LuaPostArgAmountError(func_name, 1, num_args);

  LuaCheckNilValue(func_name, L, 1);
  if (num_args >= 2)
    LuaCheckNilValue(func_name, L, 2);

  const char* temp = lua_tostring(L, 1);
  double scale = 1.0;
  if (num_args >= 2)
    scale = lua_tonumber(L, 2);
  auto new_object = new opensn::UnpartitionedMesh;

  opensn::UnpartitionedMesh::Options options;
  options.file_name = std::string(temp);
  options.scale = scale;

  new_object->ReadFromEnsightGold(options);

  opensn::unpartitionedmesh_stack.emplace_back(new_object);

  lua_pushnumber(L, static_cast<lua_Number>(opensn::unpartitionedmesh_stack.size() - 1));

  return 1;
}

int
UnpartitionedMeshFromWavefrontOBJ(lua_State* L)
{
  const std::string func_name = __FUNCTION__;
  int num_args = lua_gettop(L);
  if (num_args < 1)
    LuaPostArgAmountError(func_name, 1, num_args);

  LuaCheckNilValue(func_name, L, 1);

  const char* temp = lua_tostring(L, 1);

  auto new_object = new opensn::UnpartitionedMesh;

  opensn::UnpartitionedMesh::Options options;
  options.file_name = std::string(temp);

  new_object->ReadFromWavefrontOBJ(options);

  opensn::unpartitionedmesh_stack.emplace_back(new_object);

  lua_pushnumber(L, static_cast<lua_Number>(opensn::unpartitionedmesh_stack.size() - 1));

  return 1;
}

int
UnpartitionedMeshFromMshFormat(lua_State* L)
{
  const std::string func_name = __FUNCTION__;
  int num_args = lua_gettop(L);
  if (num_args < 1)
    LuaPostArgAmountError(func_name, 1, num_args);

  LuaCheckNilValue(func_name, L, 1);

  const char* temp = lua_tostring(L, 1);

  auto new_object = new opensn::UnpartitionedMesh;

  opensn::UnpartitionedMesh::Options options;
  options.file_name = std::string(temp);

  new_object->ReadFromMsh(options);

  opensn::unpartitionedmesh_stack.emplace_back(new_object);

  lua_pushnumber(L, static_cast<lua_Number>(opensn::unpartitionedmesh_stack.size() - 1));

  return 1;
}

int
UnpartitionedMeshFromExodusII(lua_State* L)
{
  const std::string func_name = __FUNCTION__;
  int num_args = lua_gettop(L);
  if (num_args < 1)
    LuaPostArgAmountError(func_name, 1, num_args);

  LuaCheckNilValue(func_name, L, 1);
  if (num_args >= 2)
    LuaCheckNilValue(func_name, L, 2);

  const char* temp = lua_tostring(L, 1);
  double scale = 1.0;
  if (num_args >= 2)
    scale = lua_tonumber(L, 2);
  auto new_object = new opensn::UnpartitionedMesh;

  opensn::UnpartitionedMesh::Options options;
  options.file_name = std::string(temp);
  options.scale = scale;

  new_object->ReadFromExodus(options);

  opensn::unpartitionedmesh_stack.emplace_back(new_object);

  lua_pushnumber(L, static_cast<lua_Number>(opensn::unpartitionedmesh_stack.size() - 1));

  return 1;
}

} // namespace opensnlua
