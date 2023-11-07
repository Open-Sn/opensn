#include "framework/chi_lua.h"

#include "framework/mesh/FieldFunctionInterpolation/chi_ffinterpolation.h"

#include "framework/chi_runtime.h"
#include "ffinterpol_lua.h"
#include "framework/console/chi_console.h"

RegisterLuaFunctionAsIs(chiFFInterpolationInitialize);
RegisterLuaFunctionAsIs(chiFFInterpolationExecute);

/** Initialize interpolator.
 *
\param FFIHandle int Handle to the field function interpolation.

\ingroup LuaFFInterpol
\author Jan*/
int
chiFFInterpolationInitialize(lua_State* L)
{
  const std::string fname = __FUNCTION__;
  const int num_args = lua_gettop(L);
  if (num_args != 1) LuaPostArgAmountError(fname, 1, num_args);

  // Get handle to field function
  const size_t ffihandle = lua_tonumber(L, 1);

  auto p_ffi = Chi::GetStackItemPtr(Chi::field_func_interpolation_stack, ffihandle, fname);

  p_ffi->Initialize();
  return 0;
}

/** Execute interpolator.
 *
\param FFIHandle int Handle to the field function interpolation.

\ingroup LuaFFInterpol
\author Jan*/
int
chiFFInterpolationExecute(lua_State* L)
{
  const std::string fname = __FUNCTION__;
  const int num_args = lua_gettop(L);
  if (num_args != 1) LuaPostArgAmountError(fname, 1, num_args);

  // Get handle to field function
  const size_t ffihandle = lua_tonumber(L, 1);

  auto p_ffi = Chi::GetStackItemPtr(Chi::field_func_interpolation_stack, ffihandle, fname);

  p_ffi->Execute();
  return 0;
}
