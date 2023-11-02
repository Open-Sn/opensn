#include "chi_ffinter_volume.h"

#include "chi_runtime.h"
#include "console/chi_console.h"

//###################################################################
/**Calls the designated lua function*/
double
chi_mesh::FieldFunctionInterpolationVolume::CallLuaFunction(double ff_value, int mat_id) const
{
  lua_State* L = Chi::console.GetConsoleState();
  double ret_val = 0.0;

  lua_getglobal(L, op_lua_func_.c_str());
  lua_pushnumber(L, ff_value);
  lua_pushnumber(L, mat_id);

  // 2 arguments, 1 result, 0=original error object
  if (lua_pcall(L, 2, 1, 0) == 0) { ret_val = lua_tonumber(L, -1); }
  lua_pop(L, 1);

  return ret_val;
}