#include "framework/lua.h"
#include "framework/mesh/field_function_interpolation/ffinter_point.h"
#include "framework/mesh/field_function_interpolation/ffinter_slice.h"
#include "framework/mesh/field_function_interpolation/ffinter_line.h"
#include "framework/mesh/field_function_interpolation/ffinter_volume.h"

#include "framework/runtime.h"
#include "framework/logging/log.h"

#define scint(x) static_cast<int>(x)

#include "ffinterpol_lua.h"
#include "framework/console/console.h"

using namespace opensn;

RegisterLuaFunctionAsIs(chiFFInterpolationCreate);
RegisterLuaConstantAsIs(SLICE, Varying(1));
RegisterLuaConstantAsIs(LINE, Varying(2));
RegisterLuaConstantAsIs(VOLUME, Varying(3));
RegisterLuaConstantAsIs(POINT, Varying(4));

int
chiFFInterpolationCreate(lua_State* L)
{
  auto& cur_hndlr = GetCurrentHandler();

  // Process types
  int ffitype = lua_tonumber(L, 1);
  if (ffitype == scint(FieldFunctionInterpolationType::POINT))
  {
    auto new_ffi = new FieldFunctionInterpolationPoint;

    opensn::Chi::field_func_interpolation_stack.emplace_back(new_ffi);
    const size_t index = opensn::Chi::field_func_interpolation_stack.size() - 1;
    opensn::log.LogAllVerbose2() << "Created point Field Function Interpolation";
    lua_pushnumber(L, static_cast<lua_Number>(index));
    return 1;
  }
  else if (ffitype == scint(FieldFunctionInterpolationType::SLICE))
  {
    auto new_ffi = new FieldFunctionInterpolationSlice;

    opensn::Chi::field_func_interpolation_stack.emplace_back(new_ffi);
    const size_t index = opensn::Chi::field_func_interpolation_stack.size() - 1;
    opensn::log.LogAllVerbose2() << "Created slice Field Function Interpolation";
    lua_pushnumber(L, static_cast<lua_Number>(index));
    return 1;
  }
  else if (ffitype == scint(FieldFunctionInterpolationType::LINE))
  {
    auto new_ffi = new FieldFunctionInterpolationLine;

    opensn::Chi::field_func_interpolation_stack.emplace_back(new_ffi);
    const size_t index = opensn::Chi::field_func_interpolation_stack.size() - 1;
    opensn::log.LogAllVerbose2() << "Created line Field Function Interpolation";
    lua_pushnumber(L, static_cast<lua_Number>(index));
    return 1;
  }
  else if (ffitype == scint(FieldFunctionInterpolationType::VOLUME))
  {
    auto new_ffi = new FieldFunctionInterpolationVolume;

    opensn::Chi::field_func_interpolation_stack.emplace_back(new_ffi);
    const size_t index = opensn::Chi::field_func_interpolation_stack.size() - 1;
    opensn::log.LogAllVerbose2() << "Created Volume Field Function Interpolation";
    lua_pushnumber(L, static_cast<lua_Number>(index));
    return 1;
  }
  else // Fall back
  {
    opensn::log.LogAllError() << "Invalid FFITypeIndex used in chiFFInterpolationCreate.";
    opensn::Exit(EXIT_FAILURE);
  }
  return 0;
}
