#include "framework/chi_lua.h"

#include "framework/mesh/UnpartitionedMesh/unpartitioned_mesh.h"
#include "framework/chi_runtime.h"
#include "framework/logging/chi_log.h"

#include "unpartition_mesh_lua_utils.h"
#include "framework/console/chi_console.h"

namespace chi_mesh::unpartition_mesh_lua_utils
{

RegisterLuaFunctionAsIs(chiCreateEmptyUnpartitionedMesh);
RegisterLuaFunctionAsIs(chiDestroyUnpartitionedMesh);

RegisterLuaFunctionAsIs(chiUnpartitionedMeshFromVTU);
RegisterLuaFunctionAsIs(chiUnpartitionedMeshFromPVTU);
RegisterLuaFunctionAsIs(chiUnpartitionedMeshFromEnsightGold);
RegisterLuaFunctionAsIs(chiUnpartitionedMeshFromWavefrontOBJ);
RegisterLuaFunctionAsIs(chiUnpartitionedMeshFromMshFormat);
RegisterLuaFunctionAsIs(chiUnpartitionedMeshFromExodusII);

/**Creates an empty unpartitioned mesh. An empty unpartitioned mesh
 * is meant to be manipulated with calls to chiUnpartitionedMeshUploadVertex()
 * and chiUnpartitionedMeshUploadCell(). It essentially supports building a mesh
 * manually.
 *
##_

###Example
Example usage
\code
umesh = chiCreateEmptyUnpartitionedMesh()
\endcode

\ingroup LuaUnpartitionedMesh*/
int
chiCreateEmptyUnpartitionedMesh(lua_State* L)
{
  const std::string func_name = __FUNCTION__;

  Chi::unpartitionedmesh_stack.emplace_back(new chi_mesh::UnpartitionedMesh());

  lua_pushnumber(L, static_cast<lua_Number>(Chi::unpartitionedmesh_stack.size() - 1));

  return 1;
}

/**Destroy an unpartitioned mesh. This routine should be called for
 * memory sensitive simulations because each process will have a full
 * copy of this data.
 *
\param handle int Handle to mesh.

##_

###Example
Example usage
\code
chiDestroyUnpartitionedMesh(umesh)
\endcode

\ingroup LuaUnpartitionedMesh*/
int
chiDestroyUnpartitionedMesh(lua_State* L)
{
  const std::string func_name = __FUNCTION__;
  const int num_args = lua_gettop(L);
  if (num_args != 1) LuaPostArgAmountError(func_name, 1, num_args);

  LuaCheckNilValue(func_name, L, 1);

  const int handle = lua_tointeger(L, 1);

  auto mesh_ptr = Chi::GetStackItemPtr(Chi::unpartitionedmesh_stack, handle, func_name);

  mesh_ptr->CleanUp();
  Chi::unpartitionedmesh_stack[handle] = nullptr;

  Chi::log.Log() << "Unpartitioned mesh destroyed. Memory in use = " << Chi::GetMemoryUsageInMB()
                 << " MB";
  return 0;
}

/**Creates an unpartitioned mesh from VTK Unstructured mesh files.

\param file_name char Filename of the .vtu file.
\param field char Name of the cell data field from which to read
                  material and boundary identifiers (optional).

\ingroup LuaUnpartitionedMesh

##_

### Example
An example mesh creation below:
\code
chiMeshHandlerCreate()

umesh = chiUnpartitionedMeshFromVTU("ZMeshTest_0.vtu")

chiSurfaceMesherCreate(SURFACEMESHER_PREDEFINED)
chiVolumeMesherCreate(VOLUMEMESHER_UNPARTITIONED, umesh)

chiSurfaceMesherExecute()
chiVolumeMesherExecute()
\endcode


\return Handle A handle to the newly created UnpartitionedMesh*/
int
chiUnpartitionedMeshFromVTU(lua_State* L)
{
  const std::string func_name = __FUNCTION__;
  int num_args = lua_gettop(L);
  if (num_args < 1) LuaPostArgAmountError(func_name, 1, num_args);

  LuaCheckNilValue(func_name, L, 1);
  if (num_args >= 2) LuaCheckNilValue(func_name, L, 2);

  const char* temp = lua_tostring(L, 1);
  const char* field = "";
  if (num_args >= 2) field = lua_tostring(L, 2);
  auto new_object = new chi_mesh::UnpartitionedMesh;

  chi_mesh::UnpartitionedMesh::Options options;
  options.file_name = std::string(temp);
  options.material_id_fieldname = field;
  options.boundary_id_fieldname = field;

  new_object->ReadFromVTU(options);

  Chi::unpartitionedmesh_stack.emplace_back(new_object);

  lua_pushnumber(L, static_cast<lua_Number>(Chi::unpartitionedmesh_stack.size() - 1));

  return 1;
}

/**Creates an unpartitioned mesh from VTK Partitioned Unstructured mesh files
 * (.pvtu).

\param file_name char Filename of the .vtu file.
\param field char Name of the cell data field from which to read
                  material and boundary identifiers (optional).

\ingroup LuaUnpartitionedMesh

##_

### Example
An example mesh creation below:
\code
chiMeshHandlerCreate()

umesh = chiUnpartitionedMeshFromPVTU("ZMeshTest_0.vtu")

chiSurfaceMesherCreate(SURFACEMESHER_PREDEFINED)
chiVolumeMesherCreate(VOLUMEMESHER_UNPARTITIONED, umesh)

chiSurfaceMesherExecute()
chiVolumeMesherExecute()
\endcode


\return Handle A handle to the newly created UnpartitionedMesh*/
int
chiUnpartitionedMeshFromPVTU(lua_State* L)
{
  const std::string func_name = __FUNCTION__;
  int num_args = lua_gettop(L);
  if (num_args < 1) LuaPostArgAmountError(func_name, 1, num_args);

  LuaCheckNilValue(func_name, L, 1);
  if (num_args >= 2) LuaCheckNilValue(func_name, L, 2);

  const char* temp = lua_tostring(L, 1);
  const char* field = "";
  if (num_args >= 2) field = lua_tostring(L, 2);
  auto new_object = new chi_mesh::UnpartitionedMesh;

  chi_mesh::UnpartitionedMesh::Options options;
  options.file_name = std::string(temp);
  options.material_id_fieldname = field;
  options.boundary_id_fieldname = field;

  new_object->ReadFromPVTU(options);

  Chi::unpartitionedmesh_stack.emplace_back(new_object);

  lua_pushnumber(L, static_cast<lua_Number>(Chi::unpartitionedmesh_stack.size() - 1));

  return 1;
}

/**Creates an unpartitioned mesh from starccm+ exported
Ensight Gold mesh files.

\param file_name char Filename of the .case file.
\param scale float Scale to apply to the mesh

\ingroup LuaUnpartitionedMesh

##_

### Example
An example mesh creation below:
\code
chiMeshHandlerCreate()

umesh = chiUnpartitionedMeshFromEnsightGold("resources/TestObjects/Sphere.case")

chiSurfaceMesherCreate(SURFACEMESHER_PREDEFINED)
chiVolumeMesherCreate(VOLUMEMESHER_UNPARTITIONED, umesh)

chiSurfaceMesherExecute()
chiVolumeMesherExecute()
\endcode

\return Handle A handle to the newly created UnpartitionedMesh*/
int
chiUnpartitionedMeshFromEnsightGold(lua_State* L)
{
  const std::string func_name = __FUNCTION__;
  int num_args = lua_gettop(L);
  if (num_args < 1) LuaPostArgAmountError(func_name, 1, num_args);

  LuaCheckNilValue(func_name, L, 1);
  if (num_args >= 2) LuaCheckNilValue(func_name, L, 2);

  const char* temp = lua_tostring(L, 1);
  double scale = 1.0;
  if (num_args >= 2) scale = lua_tonumber(L, 2);
  auto new_object = new chi_mesh::UnpartitionedMesh;

  chi_mesh::UnpartitionedMesh::Options options;
  options.file_name = std::string(temp);
  options.scale = scale;

  new_object->ReadFromEnsightGold(options);

  Chi::unpartitionedmesh_stack.emplace_back(new_object);

  lua_pushnumber(L, static_cast<lua_Number>(Chi::unpartitionedmesh_stack.size() - 1));

  return 1;
}

/**Creates an unpartitioned mesh from a wavefront .obj file.

\param file_name char Filename of the .case file.

\ingroup LuaUnpartitionedMesh

##_

### Example
An example mesh creation below:
\code
chiMeshHandlerCreate()

umesh =
chiUnpartitionedMeshFromWavefrontOBJ("resources/TestObjects/TriangleMesh2x2.obj")

chiSurfaceMesherCreate(SURFACEMESHER_PREDEFINED)
chiVolumeMesherCreate(VOLUMEMESHER_UNPARTITIONED, umesh)

chiSurfaceMesherExecute()
chiVolumeMesherExecute()
\endcode

\return Handle A handle to the newly created UnpartitionedMesh*/
int
chiUnpartitionedMeshFromWavefrontOBJ(lua_State* L)
{
  const std::string func_name = __FUNCTION__;
  int num_args = lua_gettop(L);
  if (num_args < 1) LuaPostArgAmountError(func_name, 1, num_args);

  LuaCheckNilValue(func_name, L, 1);

  const char* temp = lua_tostring(L, 1);

  auto new_object = new chi_mesh::UnpartitionedMesh;

  chi_mesh::UnpartitionedMesh::Options options;
  options.file_name = std::string(temp);

  new_object->ReadFromWavefrontOBJ(options);

  Chi::unpartitionedmesh_stack.emplace_back(new_object);

  lua_pushnumber(L, static_cast<lua_Number>(Chi::unpartitionedmesh_stack.size() - 1));

  return 1;
}

/**Creates an unpartitioned mesh from a .msh file.

\param file_name char Filename of the .msh file.

\ingroup LuaUnpartitionedMesh

##_

### Example
An example mesh creation below:
\code
chiMeshHandlerCreate()

umesh = chiUnpartitionedMeshFromMshFormat("File.msh")

chiSurfaceMesherCreate(SURFACEMESHER_PREDEFINED)
chiVolumeMesherCreate(VOLUMEMESHER_UNPARTITIONED, umesh)

chiSurfaceMesherExecute()
chiVolumeMesherExecute()
\endcode


\return Handle A handle to the newly created UnpartitionedMesh*/
int
chiUnpartitionedMeshFromMshFormat(lua_State* L)
{
  const std::string func_name = __FUNCTION__;
  int num_args = lua_gettop(L);
  if (num_args < 1) LuaPostArgAmountError(func_name, 1, num_args);

  LuaCheckNilValue(func_name, L, 1);

  const char* temp = lua_tostring(L, 1);

  auto new_object = new chi_mesh::UnpartitionedMesh;

  chi_mesh::UnpartitionedMesh::Options options;
  options.file_name = std::string(temp);

  new_object->ReadFromMsh(options);

  Chi::unpartitionedmesh_stack.emplace_back(new_object);

  lua_pushnumber(L, static_cast<lua_Number>(Chi::unpartitionedmesh_stack.size() - 1));

  return 1;
}

/**Creates an unpartitioned mesh from ExodusII format.

\param file_name char Filename of the .case file.
\param scale float Scale to apply to the mesh

\ingroup LuaUnpartitionedMesh

##_

### Example
An example mesh creation below:
\code
chiMeshHandlerCreate()

umesh = chiUnpartitionedMeshFromExodusII("resources/TestObjects/Mesh.e")

chiSurfaceMesherCreate(SURFACEMESHER_PREDEFINED)
chiVolumeMesherCreate(VOLUMEMESHER_UNPARTITIONED, umesh)

chiSurfaceMesherExecute()
chiVolumeMesherExecute()
\endcode

\return Handle A handle to the newly created UnpartitionedMesh*/
int
chiUnpartitionedMeshFromExodusII(lua_State* L)
{
  const std::string func_name = __FUNCTION__;
  int num_args = lua_gettop(L);
  if (num_args < 1) LuaPostArgAmountError(func_name, 1, num_args);

  LuaCheckNilValue(func_name, L, 1);
  if (num_args >= 2) LuaCheckNilValue(func_name, L, 2);

  const char* temp = lua_tostring(L, 1);
  double scale = 1.0;
  if (num_args >= 2) scale = lua_tonumber(L, 2);
  auto new_object = new chi_mesh::UnpartitionedMesh;

  chi_mesh::UnpartitionedMesh::Options options;
  options.file_name = std::string(temp);
  options.scale = scale;

  new_object->ReadFromExodus(options);

  Chi::unpartitionedmesh_stack.emplace_back(new_object);

  lua_pushnumber(L, static_cast<lua_Number>(Chi::unpartitionedmesh_stack.size() - 1));

  return 1;
}

} // namespace chi_mesh::unpartition_mesh_lua_utils
