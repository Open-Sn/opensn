#include "framework/chi_lua.h"

#include "framework/mesh/chi_mesh.h"
#include "framework/mesh/MeshHandler/chi_meshhandler.h"

#include "framework/chi_runtime.h"
#include "framework/logging/chi_log.h"

#include "lua_mesh_orthomacros.h"
#include "framework/console/chi_console.h"

RegisterLuaFunctionAsIs(chiMeshCreateUnpartitioned1DOrthoMesh);
RegisterLuaFunctionAsIs(chiMeshCreateUnpartitioned2DOrthoMesh);
RegisterLuaFunctionAsIs(chiMeshCreateUnpartitioned3DOrthoMesh);

// ###################################################################
/** Creates a 1D Mesh from an array of 1D vertices.

\param x_nodes array_float An Array of floating point numbers denoting
                           1D nodes along x-axis.

\return Two handles: unpartitioned-mesh, region

\ingroup LuaMeshMacros

##_

### Example
An example 1D mesh creation below:
\code
chiMeshHandlerCreate()
nodes={0.0,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1.0}
chiMeshCreateUnpartitioned1DOrthoMesh(nodes)
chiVolumeMesherSetProperty(PARTITION_TYPE,PARMETIS)
chiVolumeMesherExecute();
\endcode

 \author Nak*/
int
chiMeshCreateUnpartitioned1DOrthoMesh(lua_State* L)
{
  //=================================== Check argc
  const char func_name[] = "chiMeshCreateUnpartitioned1DOrthoMesh";
  int num_args = lua_gettop(L);
  if (num_args != 1) LuaPostArgAmountError(func_name, 1, num_args);

  //=================================== Check args table
  if (not lua_istable(L, 1))
  {
    Chi::log.LogAllError() << func_name << ": First argument found to not be an array.";
    Chi::Exit(EXIT_FAILURE);
  }

  //=================================== Decl vars
  int table_index = 0;
  int N = 0;
  std::vector<std::vector<double>> array(3);

  //=================================== Get first array
  table_index = 1;
  N = lua_rawlen(L, table_index);
  array[table_index - 1].resize(N);
  for (int k = 0; k < N; k++)
  {
    lua_pushnumber(L, k + 1);
    lua_gettable(L, table_index);

    array[table_index - 1][k] = lua_tonumber(L, -1);
    lua_pop(L, 1);
  }

  //=================================== Create mesh
  const size_t handle = chi_mesh::CreateUnpartitioned1DOrthoMesh(array[0]);

  //=================================== Push handles
  lua_pushnumber(L, static_cast<lua_Number>(handle));
  lua_pushnumber(L, 0);

  return 2;
}

// ###################################################################
/** Creates a 2D Orthogonal Mesh from arrays of 1D vertices.

\param x_nodes array_float An Array of floating point numbers denoting
                           1D nodes along x-axis.
\param y_nodes array_float An Array of floating point numbers denoting
                           1D nodes along y-axis.

\return Two handles: unpartitioned-mesh, region

\ingroup LuaMeshMacros

##_

### Example
An example 2D mesh creation below:
\code
chiMeshHandlerCreate()
nodesx={0.0,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1.0}
nodesy={0.0,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1.0}
chiMeshCreateUnpartitioned2DOrthoMesh(nodesx,nodesy)
chiVolumeMesherSetProperty(PARTITION_TYPE,PARMETIS)
chiVolumeMesherExecute();
\endcode

 \author Nak*/
int
chiMeshCreateUnpartitioned2DOrthoMesh(lua_State* L)
{
  //=================================== Check argc
  const char func_name[] = "chiMeshCreateUnpartitioned2DOrthoMesh";
  int num_args = lua_gettop(L);
  if (num_args != 2) LuaPostArgAmountError(func_name, 2, num_args);

  //=================================== Check args table
  if (not lua_istable(L, 1))
  {
    Chi::log.LogAllError() << func_name << ": First argument found to not be an array.";
    Chi::Exit(EXIT_FAILURE);
  }
  if (not lua_istable(L, 2))
  {
    Chi::log.LogAllError() << func_name << ": Second argument found to not be an array.";
    Chi::Exit(EXIT_FAILURE);
  }

  //=================================== Decl vars
  int table_index = 0;
  int N = 0;
  std::vector<std::vector<double>> array(3);

  //=================================== Get first array
  table_index = 1;
  N = lua_rawlen(L, table_index);
  array[table_index - 1].resize(N);
  for (int k = 0; k < N; k++)
  {
    lua_pushnumber(L, k + 1);
    lua_gettable(L, table_index);

    array[table_index - 1][k] = lua_tonumber(L, -1);
    lua_pop(L, 1);
  }
  //=================================== Get second array
  table_index = 2;
  N = lua_rawlen(L, table_index);
  array[table_index - 1].resize(N);
  for (int k = 0; k < N; k++)
  {
    lua_pushnumber(L, k + 1);
    lua_gettable(L, table_index);

    array[table_index - 1][k] = lua_tonumber(L, -1);
    lua_pop(L, 1);
  }

  //=================================== Create mesh
  const size_t handle = chi_mesh::CreateUnpartitioned2DOrthoMesh(array[0], array[1]);

  //=================================== Push handles
  lua_pushnumber(L, static_cast<lua_Number>(handle));
  lua_pushnumber(L, 0);

  return 2;
}

// ###################################################################
/** Creates a 3D Orthogonal Mesh from arrays of 1D vertices. The
 * underlying mesher is an extruder.

\param x_nodes array_float An Array of floating point numbers denoting
                           1D nodes along x-axis.
\param y_nodes array_float An Array of floating point numbers denoting
                           1D nodes along y-axis.
\param z_nodes array_float An Array of floating point numbers denoting
                           1D nodes along z-axis.

\return Two handles: unpartitioned-mesh, region

\ingroup LuaMeshMacros

##_

### Example
An example 3D mesh creation below:
\code
chiMeshHandlerCreate()
nodesx={0.0,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1.0}
nodesy={0.0,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1.0}
nodesz={0.0,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1.0}
chiMeshCreateUnpartitioned3DOrthoMesh(nodesx,nodesy,nodesz)
chiVolumeMesherSetProperty(PARTITION_TYPE,PARMETIS)
chiVolumeMesherExecute();
\endcode

 \author Nak*/
int
chiMeshCreateUnpartitioned3DOrthoMesh(lua_State* L)
{
  //=================================== Check argc
  const char func_name[] = "chiMeshCreateUnpartitioned3DOrthoMesh";
  int num_args = lua_gettop(L);
  if (num_args != 3) LuaPostArgAmountError(func_name, 3, num_args);

  //=================================== Check args table
  if (not lua_istable(L, 1))
  {
    Chi::log.LogAllError() << func_name << ": First argument found to not be an array.";
    Chi::Exit(EXIT_FAILURE);
  }
  if (not lua_istable(L, 2))
  {
    Chi::log.LogAllError() << func_name << ": Second argument found to not be an array.";
    Chi::Exit(EXIT_FAILURE);
  }
  if (not lua_istable(L, 3))
  {
    Chi::log.LogAllError() << func_name << ": Third argument found to not be an array.";
    Chi::Exit(EXIT_FAILURE);
  }

  //=================================== Decl vars
  int table_index = 0;
  int N = 0;
  std::vector<std::vector<double>> array(3);

  //=================================== Get first array
  table_index = 1;
  N = lua_rawlen(L, table_index);
  array[table_index - 1].resize(N);
  for (int k = 0; k < N; k++)
  {
    lua_pushnumber(L, k + 1);
    lua_gettable(L, table_index);

    array[table_index - 1][k] = lua_tonumber(L, -1);
    lua_pop(L, 1);
  }
  //=================================== Get second array
  table_index = 2;
  N = lua_rawlen(L, table_index);
  array[table_index - 1].resize(N);
  for (int k = 0; k < N; k++)
  {
    lua_pushnumber(L, k + 1);
    lua_gettable(L, table_index);

    array[table_index - 1][k] = lua_tonumber(L, -1);
    lua_pop(L, 1);
  }
  //=================================== Get second array
  table_index = 3;
  N = lua_rawlen(L, table_index);
  array[table_index - 1].resize(N);
  for (int k = 0; k < N; k++)
  {
    lua_pushnumber(L, k + 1);
    lua_gettable(L, table_index);

    array[table_index - 1][k] = lua_tonumber(L, -1);
    lua_pop(L, 1);
  }

  //=================================== Create mesh
  const size_t handle = chi_mesh::CreateUnpartitioned3DOrthoMesh(array[0], array[1], array[2]);

  //=================================== Push handles
  lua_pushnumber(L, static_cast<lua_Number>(handle));
  lua_pushnumber(L, 0);

  return 2;
}
