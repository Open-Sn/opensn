#include "framework/lua.h"

#include "framework/mesh/mesh.h"
#include "framework/mesh/mesh_handler/mesh_handler.h"
#include "framework/mesh/mesh_continuum/mesh_continuum.h"

#include "framework/mesh/logical_volume/logical_volume.h"

#include "framework/runtime.h"

#include "meshcutting_lua.h"
#include "framework/console/console.h"

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
