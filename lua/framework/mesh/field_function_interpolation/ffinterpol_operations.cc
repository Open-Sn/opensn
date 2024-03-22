#include "framework/lua.h"
#include "framework/field_functions/interpolation/ffinterpolation.h"
#include "framework/runtime.h"
#include "ffinterpol_lua.h"
#include "framework/console/console.h"

namespace opensnlua
{

RegisterLuaFunctionNamespace(FFInterpolationInitialize, fieldfunc, Initialize);
RegisterLuaFunctionNamespace(FFInterpolationExecute, fieldfunc, Execute);

int
FFInterpolationInitialize(lua_State* L)
{
  const std::string fname = "fieldfunc.Initialize";
  LuaCheckArgs<size_t>(L, fname);

  // Get handle to field function
  const auto ffihandle = LuaArg<size_t>(L, 1);

  auto p_ffi = opensn::GetStackItemPtr(opensn::field_func_interpolation_stack, ffihandle, fname);

  p_ffi->Initialize();
  return 0;
}

int
FFInterpolationExecute(lua_State* L)
{
  const std::string fname = "fieldfunc.Execute";
  LuaCheckArgs<size_t>(L, fname);

  // Get handle to field function
  const auto ffihandle = LuaArg<size_t>(L, 1);

  auto p_ffi = opensn::GetStackItemPtr(opensn::field_func_interpolation_stack, ffihandle, fname);

  p_ffi->Execute();
  return 0;
}

} // namespace opensnlua
