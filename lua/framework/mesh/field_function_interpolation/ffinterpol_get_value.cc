#include "framework/lua.h"
#include "framework/mesh/field_function_interpolation/ffinter_point.h"
#include "framework/mesh/field_function_interpolation/ffinter_line.h"
#include "framework/mesh/field_function_interpolation/ffinter_volume.h"

#include "framework/runtime.h"
#include "framework/logging/log.h"

#define dcastPoint(x) dynamic_cast<chi_mesh::FieldFunctionInterpolationPoint&>(x)
#define dcastLine(x) dynamic_cast<chi_mesh::FieldFunctionInterpolationLine&>(x)
#define dcastVolume(x) dynamic_cast<chi_mesh::FieldFunctionInterpolationVolume&>(x)

#include "ffinterpol_lua.h"
#include "framework/console/console.h"

RegisterLuaFunctionAsIs(chiFFInterpolationGetValue);

int
chiFFInterpolationGetValue(lua_State* L)
{
  const std::string fname = __FUNCTION__;

  int num_args = lua_gettop(L);
  if (num_args != 1) LuaPostArgAmountError("chiFFInterpolationGetValue", 1, num_args);

  // Get handle to field function
  const size_t ffihandle = lua_tonumber(L, 1);

  auto p_ffi = Chi::GetStackItemPtr(Chi::field_func_interpolation_stack, ffihandle, fname);

  if (p_ffi->Type() == chi_mesh::ff_interpolation::Type::POINT)
  {
    auto& cur_ffi_point = dcastPoint(*p_ffi);
    double value = cur_ffi_point.GetPointValue();

    lua_pushnumber(L, value);
    return 1;
  }
  else if (p_ffi->Type() == chi_mesh::ff_interpolation::Type::LINE)
  {
    auto& cur_ffi_line = dcastLine(*p_ffi);

    lua_newtable(L);

    for (int ff = 0; ff < cur_ffi_line.GetFieldFunctions().size(); ff++)
    {
      lua_pushnumber(L, ff + 1);

      lua_newtable(L);
      const auto& ff_ctx = cur_ffi_line.GetFFContexts()[ff];

      for (int p = 0; p < cur_ffi_line.GetInterpolationPoints().size(); p++)
      {
        lua_pushnumber(L, p + 1);
        lua_pushnumber(L, ff_ctx.interpolation_points_values[p]);
        lua_settable(L, -3);
      }

      lua_settable(L, -3);
    }

    return 1;
  }
  else if (p_ffi->Type() == chi_mesh::ff_interpolation::Type::VOLUME)
  {
    auto& cur_ffi_volume = dcastVolume(*p_ffi);
    double value = cur_ffi_volume.GetOpValue();

    lua_pushnumber(L, value);
    return 1;
  }
  else
  {
    Chi::log.Log0Warning() << "chiFFInterpolationGetValue is currently only supported for "
                           << " POINT, LINE and VOLUME interpolator types.";
  }

  return 0;
}
