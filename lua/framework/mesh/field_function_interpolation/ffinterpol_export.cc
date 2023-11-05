#include "opensn/framework/chi_lua.h"

#include "opensn/framework/mesh/FieldFunctionInterpolation/chi_ffinterpolation.h"

#include "opensn/framework/chi_runtime.h"
#include "opensn/framework/logging/chi_log.h"
#include "ffinterpol_lua.h"
#include "opensn/framework/console/chi_console.h"

RegisterLuaFunctionAsIs(chiFFInterpolationExportPython);

//###################################################################
/** Export interpolation to python line,contour plot depending on the
 * type of interpolation.
 *
\param FFIHandle int Handle to the field function interpolation.
\param BaseName char Base name to be used for exported files.

\ingroup LuaFFInterpol
\author Jan*/
int
chiFFInterpolationExportPython(lua_State* L)
{
  const std::string fname = __FUNCTION__;

  const int num_args = lua_gettop(L);
  if (num_args < 1) LuaPostArgAmountError(fname, 1, num_args);

  //================================================== Get handle to field
  // function
  const size_t ffihandle = lua_tonumber(L, 1);

  auto p_ffi = Chi::GetStackItemPtr(Chi::field_func_interpolation_stack, ffihandle, fname);

  std::string base_name = p_ffi->GetDefaultFileBaseName() + std::to_string(ffihandle);
  if (num_args == 2) base_name = lua_tostring(L, 2);

  p_ffi->ExportPython(base_name);

  return 0;
}
