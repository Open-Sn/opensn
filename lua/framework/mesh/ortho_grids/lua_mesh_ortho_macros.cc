#include "framework/lua.h"
#include "framework/data_types/varying.h"
#include "framework/mesh/mesh_continuum/mesh_continuum.h"
#include "framework/mesh/volume_mesher/volume_mesher.h"
#include "framework/runtime.h"
#include "framework/logging/log.h"
#include "framework/utils/timer.h"
#include "lua_mesh_ortho_macros.h"
#include "framework/console/console.h"

using namespace opensn;

RegisterLuaConstant(OrthoBoundaryID, XMAX, Varying(0));
RegisterLuaConstant(OrthoBoundaryID, XMIN, Varying(1));
RegisterLuaConstant(OrthoBoundaryID, YMAX, Varying(2));
RegisterLuaConstant(OrthoBoundaryID, YMIN, Varying(3));
RegisterLuaConstant(OrthoBoundaryID, ZMAX, Varying(4));
RegisterLuaConstant(OrthoBoundaryID, ZMIN, Varying(5));

RegisterLuaFunctionAsIs(MeshCreateUnpartitioned1DOrthoMesh);
RegisterLuaFunctionAsIs(MeshCreateUnpartitioned2DOrthoMesh);
RegisterLuaFunctionAsIs(MeshCreateUnpartitioned3DOrthoMesh);
RegisterLuaFunctionNamespace(MeshSetupOrthogonalBoundaries, mesh, SetupOrthogonalBoundaries);

int
MeshCreateUnpartitioned1DOrthoMesh(lua_State* L)
{
  // Check argc
  const char func_name[] = "MeshCreateUnpartitioned1DOrthoMesh";
  int num_args = lua_gettop(L);
  if (num_args != 1)
    LuaPostArgAmountError(func_name, 1, num_args);

  // Check args table
  if (not lua_istable(L, 1))
  {
    opensn::log.LogAllError() << func_name << ": First argument found to not be an array.";
    opensn::Exit(EXIT_FAILURE);
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
MeshCreateUnpartitioned2DOrthoMesh(lua_State* L)
{
  // Check argc
  const char func_name[] = "MeshCreateUnpartitioned2DOrthoMesh";
  int num_args = lua_gettop(L);
  if (num_args != 2)
    LuaPostArgAmountError(func_name, 2, num_args);

  // Check args table
  if (not lua_istable(L, 1))
  {
    opensn::log.LogAllError() << func_name << ": First argument found to not be an array.";
    opensn::Exit(EXIT_FAILURE);
  }
  if (not lua_istable(L, 2))
  {
    opensn::log.LogAllError() << func_name << ": Second argument found to not be an array.";
    opensn::Exit(EXIT_FAILURE);
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
MeshCreateUnpartitioned3DOrthoMesh(lua_State* L)
{
  // Check argc
  const char func_name[] = "MeshCreateUnpartitioned3DOrthoMesh";
  int num_args = lua_gettop(L);
  if (num_args != 3)
    LuaPostArgAmountError(func_name, 3, num_args);

  // Check args table
  if (not lua_istable(L, 1))
  {
    opensn::log.LogAllError() << func_name << ": First argument found to not be an array.";
    opensn::Exit(EXIT_FAILURE);
  }
  if (not lua_istable(L, 2))
  {
    opensn::log.LogAllError() << func_name << ": Second argument found to not be an array.";
    opensn::Exit(EXIT_FAILURE);
  }
  if (not lua_istable(L, 3))
  {
    opensn::log.LogAllError() << func_name << ": Third argument found to not be an array.";
    opensn::Exit(EXIT_FAILURE);
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

int
MeshSetupOrthogonalBoundaries(lua_State* L)
{
  opensn::log.Log() << program_timer.GetTimeString() << " Setting orthogonal boundaries.";

  auto vol_cont = GetCurrentMesh();

  const Vector3 ihat(1.0, 0.0, 0.0);
  const Vector3 jhat(0.0, 1.0, 0.0);
  const Vector3 khat(0.0, 0.0, 1.0);

  for (auto& cell : vol_cont->local_cells)
    for (auto& face : cell.faces_)
      if (not face.has_neighbor_)
      {
        Vector3& n = face.normal_;

        std::string boundary_name;
        if (n.Dot(ihat) > 0.999)
          boundary_name = "XMAX";
        else if (n.Dot(ihat) < -0.999)
          boundary_name = "XMIN";
        else if (n.Dot(jhat) > 0.999)
          boundary_name = "YMAX";
        else if (n.Dot(jhat) < -0.999)
          boundary_name = "YMIN";
        else if (n.Dot(khat) > 0.999)
          boundary_name = "ZMAX";
        else if (n.Dot(khat) < -0.999)
          boundary_name = "ZMIN";

        uint64_t bndry_id = vol_cont->MakeBoundaryID(boundary_name);

        face.neighbor_id_ = bndry_id;

        vol_cont->GetBoundaryIDMap()[bndry_id] = boundary_name;
      } // if bndry

  opensn::mpi_comm.barrier();
  opensn::log.Log() << program_timer.GetTimeString() << " Done setting orthogonal boundaries.";

  return 0;
}
