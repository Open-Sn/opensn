#include "framework/lua.h"
#include "unpartition_mesh_lua_utils.h"

#include "framework/runtime.h"
#include "framework/logging/log.h"
#include "framework/console/console.h"

namespace chi_mesh::unpartition_mesh_lua_utils
{

RegisterLuaFunctionAsIs(chiUnpartitionedMeshUploadVertex);
RegisterLuaFunctionAsIs(chiUnpartitionedMeshUploadCell);
RegisterLuaFunctionAsIs(chiUnpartitionedMeshFinalizeEmpty);

int
chiUnpartitionedMeshUploadVertex(lua_State* L)
{
  const std::string fname = __FUNCTION__;
  const int num_args = lua_gettop(L);
  if (num_args != 4) LuaPostArgAmountError(fname, 4, num_args);

  LuaCheckNilValue(fname, L, 1);
  LuaCheckNilValue(fname, L, 2);
  LuaCheckNilValue(fname, L, 3);
  LuaCheckNilValue(fname, L, 4);

  const int handle = lua_tointeger(L, 1);
  const double x = lua_tonumber(L, 2);
  const double y = lua_tonumber(L, 3);
  const double z = lua_tonumber(L, 4);

  auto& mesh =
    Chi::GetStackItem<chi_mesh::UnpartitionedMesh>(Chi::unpartitionedmesh_stack, handle, fname);

  mesh.GetVertices().emplace_back(x, y, z);

  size_t vert_handle = mesh.GetVertices().size() - 1;

  lua_pushinteger(L, static_cast<lua_Integer>(vert_handle));
  return 1;
}

int
chiUnpartitionedMeshUploadCell(lua_State* L)
{
  const std::string fname = __FUNCTION__;
  const int num_args = lua_gettop(L);
  if (num_args < 2) LuaPostArgAmountError(fname, 2, num_args);

  LuaCheckNilValue(fname, L, 1);
  LuaCheckNilValue(fname, L, 2);
  if (num_args == 3) LuaCheckBoolValue(fname, L, 3);

  const int handle = lua_tointeger(L, 1);
  bool verbose = false;
  if (num_args == 3) verbose = lua_toboolean(L, 3);

  auto& mesh =
    Chi::GetStackItem<chi_mesh::UnpartitionedMesh>(Chi::unpartitionedmesh_stack, handle, fname);

  LuaCheckTableValue(fname, L, 2);

  auto GetField = [fname](lua_State* L, int index, const std::string& field_name)
  {
    if (not lua_getfield(L, index, field_name.c_str()))
    {
      std::stringstream message;

      message << fname << ": Cell table could not be used because the "
              << "field \"" << field_name << "\" is missing.";

      throw std::logic_error(message.str());
    }
  };

  GetField(L, 2, "type");
  const std::string cell_type_str = lua_tostring(L, -1);
  lua_pop(L, 1);

  GetField(L, 2, "sub_type");
  const std::string cell_sub_type_str = lua_tostring(L, -1);
  lua_pop(L, 1);

  GetField(L, 2, "num_faces");
  const int cell_num_faces = lua_tointeger(L, -1);
  lua_pop(L, 1);

  int cell_material_id = -1;
  if (lua_getfield(L, 2, "material_id"))
  {
    cell_material_id = lua_tointeger(L, -1);
    lua_pop(L, 1);
  }

  if (verbose)
  {
    Chi::log.Log() << "Cell type       : " << cell_type_str;
    Chi::log.Log() << "Cell sub-type   : " << cell_sub_type_str;
    Chi::log.Log() << "Cell num_faces  : " << cell_num_faces;
    Chi::log.Log() << "Cell material_id: " << cell_material_id;
  }

  std::vector<std::vector<uint64_t>> proxy_faces(cell_num_faces);
  const std::string face_prefix = "face";
  for (int f = 0; f < cell_num_faces; ++f)
  {
    auto face_field_name = face_prefix + std::to_string(f);

    GetField(L, 2, face_field_name);

    std::vector<double> vals;
    const int table_index = lua_gettop(L);
    LuaPopulateVectorFrom1DArray(fname, L, table_index, vals);

    if (verbose)
    {
      std::stringstream outstr;
      outstr << "face" << f << " ";
      for (auto val : vals)
        outstr << val << " ";
      Chi::log.Log() << outstr.str();
    }

    std::vector<uint64_t> proxy_face;
    proxy_face.reserve(vals.size());
    for (auto val : vals)
      proxy_face.push_back(static_cast<uint64_t>(val));

    proxy_faces[f] = std::move(proxy_face);

    lua_pop(L, 1); // pop off the table
  }

  mesh.PushProxyCell(
    cell_type_str, cell_sub_type_str, cell_num_faces, cell_material_id, proxy_faces);

  size_t cell_handle = mesh.GetNumberOfCells() - 1;

  lua_pushinteger(L, static_cast<lua_Integer>(cell_handle));
  return 1;
}

int
chiUnpartitionedMeshFinalizeEmpty(lua_State* L)
{
  const std::string fname = __FUNCTION__;
  const int num_args = lua_gettop(L);
  if (num_args != 1) LuaPostArgAmountError(fname, 1, num_args);

  LuaCheckNilValue(fname, L, 1);

  const int handle = lua_tointeger(L, 1);

  auto& mesh =
    Chi::GetStackItem<chi_mesh::UnpartitionedMesh>(Chi::unpartitionedmesh_stack, handle, fname);

  mesh.ComputeCentroidsAndCheckQuality();
  mesh.BuildMeshConnectivity();

  return 0;
}

} // namespace chi_mesh::unpartition_mesh_lua_utils
