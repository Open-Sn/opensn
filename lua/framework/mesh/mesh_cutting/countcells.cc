#include "framework/chi_lua.h"

#include "framework/mesh/chi_mesh.h"
#include "framework/mesh/MeshHandler/chi_meshhandler.h"
#include "framework/mesh/MeshContinuum/chi_meshcontinuum.h"

#include "framework/mesh/LogicalVolume/LogicalVolume.h"

#include "framework/chi_runtime.h"

#include "meshcutting_lua.h"
#include "framework/console/chi_console.h"

RegisterLuaFunctionAsIs(chiCountMeshInLogicalVolume);

int
chiCountMeshInLogicalVolume(lua_State* L)
{
  const std::string fname = __FUNCTION__;

  // Arg checking
  int num_args = lua_gettop(L);
  if (num_args != 1) LuaPostArgAmountError(__FUNCTION__, 1, num_args);

  LuaCheckNilValue(__FUNCTION__, L, 1);

  int log_vol_handle = lua_tonumber(L, 1);

  auto& handler = chi_mesh::GetCurrentHandler();

  const auto& log_vol =
    Chi::GetStackItem<chi_mesh::LogicalVolume>(Chi::object_stack, log_vol_handle, fname);

  auto& grid = handler.GetGrid();

  size_t count = grid->CountCellsInLogicalVolume(log_vol);

  lua_pushinteger(L, int(count));
  return 1;
}
