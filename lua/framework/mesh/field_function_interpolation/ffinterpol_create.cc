#include "framework/chi_lua.h"
#include "framework/mesh/FieldFunctionInterpolation/Point/chi_ffinter_point.h"
#include "framework/mesh/FieldFunctionInterpolation/Slice/chi_ffinter_slice.h"
#include "framework/mesh/FieldFunctionInterpolation/Line/chi_ffinter_line.h"
#include "framework/mesh/FieldFunctionInterpolation/Volume/chi_ffinter_volume.h"

#include "framework/chi_runtime.h"
#include "framework/logging/chi_log.h"

#define scint(x) static_cast<int>(x)

#include "ffinterpol_lua.h"
#include "framework/console/chi_console.h"

RegisterLuaFunctionAsIs(chiFFInterpolationCreate);
RegisterLuaConstantAsIs(SLICE, chi_data_types::Varying(1));
RegisterLuaConstantAsIs(LINE, chi_data_types::Varying(2));
RegisterLuaConstantAsIs(VOLUME, chi_data_types::Varying(3));
RegisterLuaConstantAsIs(POINT, chi_data_types::Varying(4));

// #############################################################################
/** Creates a new field function interpolation.
 *
\param FFITypeIndex int Type of field function interpolation.

##_

###FFITypeIndex\n
POINT           = A point probe. \n
SLICE           = Two dimensional slice of the mesh. \n
LINE            = Line defined by two points.\n
VOLUME          = Volume either referring to the entire volume or that of a
                  logical volume assigned to the interpolator.\n

\return Handle int Handle to the created interpolation.
\ingroup LuaFFInterpol
\author Jan*/
int
chiFFInterpolationCreate(lua_State* L)
{
  auto& cur_hndlr = chi_mesh::GetCurrentHandler();

  //================================================== Process types
  int ffitype = lua_tonumber(L, 1);
  if (ffitype == scint(chi_mesh::ff_interpolation::Type::POINT))
  {
    auto new_ffi = new chi_mesh::FieldFunctionInterpolationPoint;

    Chi::field_func_interpolation_stack.emplace_back(new_ffi);
    const size_t index = Chi::field_func_interpolation_stack.size() - 1;
    Chi::log.LogAllVerbose2() << "Created point Field Function Interpolation";
    lua_pushnumber(L, static_cast<lua_Number>(index));
    return 1;
  }
  else if (ffitype == scint(chi_mesh::ff_interpolation::Type::SLICE))
  {
    auto new_ffi = new chi_mesh::FieldFunctionInterpolationSlice;

    Chi::field_func_interpolation_stack.emplace_back(new_ffi);
    const size_t index = Chi::field_func_interpolation_stack.size() - 1;
    Chi::log.LogAllVerbose2() << "Created slice Field Function Interpolation";
    lua_pushnumber(L, static_cast<lua_Number>(index));
    return 1;
  }
  else if (ffitype == scint(chi_mesh::ff_interpolation::Type::LINE))
  {
    auto new_ffi = new chi_mesh::FieldFunctionInterpolationLine;

    Chi::field_func_interpolation_stack.emplace_back(new_ffi);
    const size_t index = Chi::field_func_interpolation_stack.size() - 1;
    Chi::log.LogAllVerbose2() << "Created line Field Function Interpolation";
    lua_pushnumber(L, static_cast<lua_Number>(index));
    return 1;
  }
  else if (ffitype == scint(chi_mesh::ff_interpolation::Type::VOLUME))
  {
    auto new_ffi = new chi_mesh::FieldFunctionInterpolationVolume;

    Chi::field_func_interpolation_stack.emplace_back(new_ffi);
    const size_t index = Chi::field_func_interpolation_stack.size() - 1;
    Chi::log.LogAllVerbose2() << "Created Volume Field Function Interpolation";
    lua_pushnumber(L, static_cast<lua_Number>(index));
    return 1;
  }
  else // Fall back
  {
    Chi::log.LogAllError() << "Invalid FFITypeIndex used in chiFFInterpolationCreate.";
    Chi::Exit(EXIT_FAILURE);
  }
  return 0;
}
