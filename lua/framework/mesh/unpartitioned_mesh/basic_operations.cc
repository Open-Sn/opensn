#include "opensn/framework/chi_lua.h"
#include "unpartition_mesh_lua_utils.h"

#include "opensn/framework/chi_runtime.h"
#include "opensn/framework/logging/chi_log.h"
#include "opensn/framework/console/chi_console.h"

namespace chi_mesh::unpartition_mesh_lua_utils
{

RegisterLuaFunctionAsIs(chiUnpartitionedMeshUploadVertex);
RegisterLuaFunctionAsIs(chiUnpartitionedMeshUploadCell);
RegisterLuaFunctionAsIs(chiUnpartitionedMeshFinalizeEmpty);

//###################################################################
/**Uploads a vertex.
 *
\param handle int Handle to mesh.
\param x      double x-coordinate.
\param y      double y-coordinate.
\param z      double z-coordinate.

## _

###Example
Example usage
\code
chiUnpartitionedMeshUploadVertex(umesh, 0, 0, 0)
chiUnpartitionedMeshUploadVertex(umesh, 1, 0, 0)
chiUnpartitionedMeshUploadVertex(umesh, 1, 1, 0)
chiUnpartitionedMeshUploadVertex(umesh, 0, 1, 0)
\endcode

\ingroup LuaUnpartitionedMesh
 */
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

//###################################################################
/**Uploads a cell

\param handle int Handle to mesh.
\param cell_table lua_table A Lua-table containing fields of data. See
 cell_table below.

## _

###cell_table
A lua-table with the following fields:
- `type` <I>string</I>. The value of this field contains the cell's primary
type. The value can be "SLAB", "POLYGON" or "POLYHEDRON".
- `sub_type` <I>string</I>. The value of this field constains the cell's
 secondary type. The value can be "SLAB", "POLYGON", "TRIANGLE",
 "QUADRILATERAL", "POLYHEDRON", "TETRAHEDRON" or "HEXAHEDRON".
- `num_faces` <I>int</I>. The value of this field represent the number of faces
 specified for this cell. Each face is contained in a field "faceX" where X is
 the face index.
- `material_id` <I>int</I>. (Optional) The value of this field holds a material
 identifier. If not provided, will be defaulted to -1.
- `faceX` <I>table</I>. A field holding a lua-table containing the vertex-ids
 of face X. There must be `num_faces` of these fields, i.e., "face0", "face1",
etc.

###Example
Example usage
\code
chiUnpartitionedMeshUploadVertex(umesh, 0, 0, 0)
chiUnpartitionedMeshUploadVertex(umesh, 1, 0, 0)
chiUnpartitionedMeshUploadVertex(umesh, 1, 1, 0)
chiUnpartitionedMeshUploadVertex(umesh, 0, 1, 0)

chiUnpartitionedMeshUploadVertex(umesh, 0, 0, 1)
chiUnpartitionedMeshUploadVertex(umesh, 1, 0, 1)
chiUnpartitionedMeshUploadVertex(umesh, 1, 1, 1)
chiUnpartitionedMeshUploadVertex(umesh, 0, 1, 1)

cell = {}
cell.type        = "POLYHEDRON"
cell.sub_type    = "HEXAHEDRON"
cell.num_faces   = 6
cell.material_id = 0
cell.face0 = {1,2,6,5}
cell.face1 = {0,4,7,3}
cell.face2 = {2,3,7,6}
cell.face3 = {0,1,5,4}
cell.face4 = {4,5,6,7}
cell.face5 = {0,3,2,1}

chiUnpartitionedMeshUploadCell(umesh, cell, true)
chiUnpartitionedMeshFinalizeEmpty(umesh)
\endcode

\ingroup LuaUnpartitionedMesh
 */
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

//###################################################################
/**Finalizes a mesh. This usually involves computing centroids and
 * establishing connectivity.
 *
\param handle int Handle to mesh.

## _

###Example
Example usage
\code
chiUnpartitionedMeshUploadVertex(umesh, 0, 0, 0)
chiUnpartitionedMeshUploadVertex(umesh, 1, 0, 0)
chiUnpartitionedMeshUploadVertex(umesh, 1, 1, 0)
chiUnpartitionedMeshUploadVertex(umesh, 0, 1, 0)

chiUnpartitionedMeshUploadVertex(umesh, 0, 0, 1)
chiUnpartitionedMeshUploadVertex(umesh, 1, 0, 1)
chiUnpartitionedMeshUploadVertex(umesh, 1, 1, 1)
chiUnpartitionedMeshUploadVertex(umesh, 0, 1, 1)

cell = {}
cell.type        = "POLYHEDRON"
cell.sub_type    = "HEXAHEDRON"
cell.num_faces   = 6
cell.material_id = 0
cell.face0 = {1,2,6,5}
cell.face1 = {0,4,7,3}
cell.face2 = {2,3,7,6}
cell.face3 = {0,1,5,4}
cell.face4 = {4,5,6,7}
cell.face5 = {0,3,2,1}

chiUnpartitionedMeshUploadCell(umesh, cell, true)
chiUnpartitionedMeshFinalizeEmpty(umesh)
\endcode

\ingroup LuaUnpartitionedMesh
 */
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
