// SPDX-FileCopyrightText: 2024 The OpenSn Authors <https://open-sn.github.io/opensn/>
// SPDX-License-Identifier: MIT

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
  const std::string fname = "fieldfunc.ExportPython";
  LuaCheckArgs<int>(L, fname);

  // Get handle to field function
  const auto ffihandle = LuaArg<size_t>(L, 1);

  auto p_ffi = opensn::GetStackItemPtr(opensn::field_func_interpolation_stack, ffihandle, fname);

  auto base_name =
    LuaArgOptional<std::string>(L, 2, p_ffi->GetDefaultFileBaseName() + std::to_string(ffihandle));

  p_ffi->ExportPython(base_name);

  return LuaReturn(L);
}

} // namespace opensnlua
