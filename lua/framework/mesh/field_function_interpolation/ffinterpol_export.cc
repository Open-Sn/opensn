// SPDX-FileCopyrightText: 2024 The OpenSn Authors <https://open-sn.github.io/opensn/>
// SPDX-License-Identifier: MIT

#include "framework/lua.h"
#include "ffinterpol_lua.h"
#include "framework/field_functions/interpolation/ffinterpolation.h"
#include "framework/console/console.h"
#include "framework/logging/log.h"
#include "framework/runtime.h"

namespace opensnlua
{

RegisterLuaFunctionInNamespace(FFInterpolationExport, fieldfunc, Export);

int
FFInterpolationExport(lua_State* L)
{
  LuaCheckArgs<int>(L, "fieldfunc.Export");

  // Get handle to field function
  const auto ffihandle = LuaArg<size_t>(L, 1);
  auto p_ffi =
    opensn::GetStackItemPtr(opensn::field_func_interpolation_stack, ffihandle, "fieldfunc.Export");
  auto base_name = LuaArgOptional<std::string>(L, 2, opensn::input_path.stem());
  p_ffi->Export(base_name);

  return LuaReturn(L);
}

} // namespace opensnlua
