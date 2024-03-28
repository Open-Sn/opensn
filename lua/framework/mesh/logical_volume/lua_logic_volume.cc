#include "lua_logic_volume.h"
#include "framework/mesh/logical_volume/logical_volume.h"
#include "framework/runtime.h"
#include "framework/console/console.h"

using namespace opensn;

namespace opensnlua
{

RegisterLuaFunctionNamespace(LogVolPointSense, logvol, PointSense);

int
LogVolPointSense(lua_State* L)
{
  const std::string fname = "LogVolPointSense";
  const int num_args = lua_gettop(L);
  if (num_args != 4)
    LuaPostArgAmountError(fname, 4, num_args);

  LuaCheckNilValue(fname, L, 1);

  const int lv_handle = lua_tointeger(L, 1);

  const auto& lv = opensn::GetStackItem<LogicalVolume>(opensn::object_stack, lv_handle, fname);

  const Vector3 point(lua_tonumber(L, 2), lua_tonumber(L, 3), lua_tonumber(L, 4));

  if (lv.Inside(point))
    lua_pushboolean(L, true);
  else
    lua_pushboolean(L, false);

  return 1;
}

} // namespace opensnlua
