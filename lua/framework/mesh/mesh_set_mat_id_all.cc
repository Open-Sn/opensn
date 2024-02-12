#include "framework/lua.h"
#include "framework/runtime.h"
#include "framework/mesh/mesh_continuum/mesh_continuum.h"
#include "framework/logging/log.h"
#include "framework/utils/timer.h"
#include "framework/console/console.h"

int MeshSetMatIDToAll(lua_State* L);

RegisterLuaFunctionNamespace(MeshSetMatIDToAll, mesh, SetMatIDToAll);

int
MeshSetMatIDToAll(lua_State* L)
{
  int num_args = lua_gettop(L);
  if (num_args != 1)
    LuaPostArgAmountError(__FUNCTION__, 1, num_args);

  LuaCheckNilValue(__FUNCTION__, L, 1);

  int mat_id = lua_tonumber(L, 1);

  auto vol_cont = opensn::GetCurrentMesh();
  vol_cont->SetMatIDToAll(mat_id);
  opensn::mpi_comm.barrier();
  opensn::log.Log() << opensn::program_timer.GetTimeString() << " Done setting material id "
                    << mat_id << " to all cells";

  return 0;
}
