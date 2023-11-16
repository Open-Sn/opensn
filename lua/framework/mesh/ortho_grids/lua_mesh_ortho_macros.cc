#include "framework/lua.h"

#include "framework/mesh/mesh.h"
#include "framework/mesh/mesh_handler/mesh_handler.h"

#include "framework/runtime.h"
#include "framework/logging/log.h"

#include "lua_mesh_ortho_macros.h"
#include "framework/console/console.h"

RegisterLuaFunctionAsIs(chiMeshCreateUnpartitioned1DOrthoMesh);
RegisterLuaFunctionAsIs(chiMeshCreateUnpartitioned2DOrthoMesh);
RegisterLuaFunctionAsIs(chiMeshCreateUnpartitioned3DOrthoMesh);

int
chiMeshCreateUnpartitioned1DOrthoMesh(lua_State* L)
{
  // Check argc
  const char func_name[] = "chiMeshCreateUnpartitioned1DOrthoMesh";
  int num_args = lua_gettop(L);
  if (num_args != 1) LuaPostArgAmountError(func_name, 1, num_args);

  // Check args table
  if (not lua_istable(L, 1))
  {
    opensn::log.LogAllError() << func_name << ": First argument found to not be an array.";
    opensn::Chi::Exit(EXIT_FAILURE);
  }

  // Decl vars
  int table_index = 0;
  int N = 0;
  std::vector<std::vector<double>> array(3);

  // Get first array
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

  // Create mesh
  const size_t handle = opensn::CreateUnpartitioned1DOrthoMesh(array[0]);

  // Push handles
  lua_pushnumber(L, static_cast<lua_Number>(handle));
  lua_pushnumber(L, 0);

  return 2;
}

int
chiMeshCreateUnpartitioned2DOrthoMesh(lua_State* L)
{
  // Check argc
  const char func_name[] = "chiMeshCreateUnpartitioned2DOrthoMesh";
  int num_args = lua_gettop(L);
  if (num_args != 2) LuaPostArgAmountError(func_name, 2, num_args);

  // Check args table
  if (not lua_istable(L, 1))
  {
    opensn::log.LogAllError() << func_name << ": First argument found to not be an array.";
    opensn::Chi::Exit(EXIT_FAILURE);
  }
  if (not lua_istable(L, 2))
  {
    opensn::log.LogAllError() << func_name << ": Second argument found to not be an array.";
    opensn::Chi::Exit(EXIT_FAILURE);
  }

  // Decl vars
  int table_index = 0;
  int N = 0;
  std::vector<std::vector<double>> array(3);

  // Get first array
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
  // Get second array
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

  // Create mesh
  const size_t handle = opensn::CreateUnpartitioned2DOrthoMesh(array[0], array[1]);

  // Push handles
  lua_pushnumber(L, static_cast<lua_Number>(handle));
  lua_pushnumber(L, 0);

  return 2;
}

int
chiMeshCreateUnpartitioned3DOrthoMesh(lua_State* L)
{
  // Check argc
  const char func_name[] = "chiMeshCreateUnpartitioned3DOrthoMesh";
  int num_args = lua_gettop(L);
  if (num_args != 3) LuaPostArgAmountError(func_name, 3, num_args);

  // Check args table
  if (not lua_istable(L, 1))
  {
    opensn::log.LogAllError() << func_name << ": First argument found to not be an array.";
    opensn::Chi::Exit(EXIT_FAILURE);
  }
  if (not lua_istable(L, 2))
  {
    opensn::log.LogAllError() << func_name << ": Second argument found to not be an array.";
    opensn::Chi::Exit(EXIT_FAILURE);
  }
  if (not lua_istable(L, 3))
  {
    opensn::log.LogAllError() << func_name << ": Third argument found to not be an array.";
    opensn::Chi::Exit(EXIT_FAILURE);
  }

  // Decl vars
  int table_index = 0;
  int N = 0;
  std::vector<std::vector<double>> array(3);

  // Get first array
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
  // Get second array
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
  // Get second array
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

  // Create mesh
  const size_t handle = opensn::CreateUnpartitioned3DOrthoMesh(array[0], array[1], array[2]);

  // Push handles
  lua_pushnumber(L, static_cast<lua_Number>(handle));
  lua_pushnumber(L, 0);

  return 2;
}
