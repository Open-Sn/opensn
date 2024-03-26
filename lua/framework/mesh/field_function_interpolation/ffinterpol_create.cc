#include "framework/lua.h"
#include "framework/field_functions/interpolation/ffinter_point.h"
#include "framework/field_functions/interpolation/ffinter_slice.h"
#include "framework/field_functions/interpolation/ffinter_line.h"
#include "framework/field_functions/interpolation/ffinter_volume.h"
#include "framework/runtime.h"
#include "framework/logging/log.h"
#include "ffinterpol_lua.h"
#include "framework/console/console.h"

using namespace opensn;

namespace opensnlua
{

RegisterLuaFunctionNamespace(FFInterpolationCreate, fieldfunc, FFInterpolationCreate);
RegisterLuaConstantAsIs(SLICE, Varying(1));
RegisterLuaConstantAsIs(LINE, Varying(2));
RegisterLuaConstantAsIs(VOLUME, Varying(3));
RegisterLuaConstantAsIs(POINT, Varying(4));

int
FFInterpolationCreate(lua_State* L)
{
  // Process types
  int ffitype = lua_tonumber(L, 1);
  if (ffitype == static_cast<int>(FieldFunctionInterpolationType::POINT))
  {
    auto new_ffi = new FieldFunctionInterpolationPoint;

    opensn::field_func_interpolation_stack.emplace_back(new_ffi);
    const size_t index = opensn::field_func_interpolation_stack.size() - 1;
    opensn::log.LogAllVerbose2() << "Created point Field Function Interpolation";
    lua_pushnumber(L, static_cast<lua_Number>(index));
    return 1;
  }
  else if (ffitype == static_cast<int>(FieldFunctionInterpolationType::SLICE))
  {
    auto new_ffi = new FieldFunctionInterpolationSlice;

    opensn::field_func_interpolation_stack.emplace_back(new_ffi);
    const size_t index = opensn::field_func_interpolation_stack.size() - 1;
    opensn::log.LogAllVerbose2() << "Created slice Field Function Interpolation";
    lua_pushnumber(L, static_cast<lua_Number>(index));
    return 1;
  }
  else if (ffitype == static_cast<int>(FieldFunctionInterpolationType::LINE))
  {
    auto new_ffi = new FieldFunctionInterpolationLine;

    opensn::field_func_interpolation_stack.emplace_back(new_ffi);
    const size_t index = opensn::field_func_interpolation_stack.size() - 1;
    opensn::log.LogAllVerbose2() << "Created line Field Function Interpolation";
    lua_pushnumber(L, static_cast<lua_Number>(index));
    return 1;
  }
  else if (ffitype == static_cast<int>(FieldFunctionInterpolationType::VOLUME))
  {
    auto new_ffi = new FieldFunctionInterpolationVolume;

    opensn::field_func_interpolation_stack.emplace_back(new_ffi);
    const size_t index = opensn::field_func_interpolation_stack.size() - 1;
    opensn::log.LogAllVerbose2() << "Created Volume Field Function Interpolation";
    lua_pushnumber(L, static_cast<lua_Number>(index));
    return 1;
  }
  else // Fall back
  {
    opensn::log.LogAllError() << "Invalid FFITypeIndex used in FFInterpolationCreate.";
    opensn::Exit(EXIT_FAILURE);
  }
  return 0;
}

} // namespace opensnlua
