#include "framework/lua.h"
#include "framework/field_functions/interpolation/ffinterpolation.h"
#include "framework/runtime.h"
#include "framework/logging/log.h"
#include "ffinterpol_lua.h"
#include "framework/console/console.h"

namespace opensnlua
{

RegisterLuaFunctionNamespace(FFInterpolationExportPython, fieldfunc, ExportPython);

int
FFInterpolationExportPython(lua_State* L)
{
  const std::string fname = __FUNCTION__;

  const int num_args = lua_gettop(L);
  if (num_args < 1)
    LuaPostArgAmountError(fname, 1, num_args);

  // Get handle to field function
  const auto ffihandle = LuaArg<size_t>(L, 1);

  auto p_ffi = opensn::GetStackItemPtr(opensn::field_func_interpolation_stack, ffihandle, fname);

  auto base_name =
    LuaArgOptional<std::string>(L, 2, p_ffi->GetDefaultFileBaseName() + std::to_string(ffihandle));

  p_ffi->ExportPython(base_name);

  return 0;
}

} // namespace opensnlua
