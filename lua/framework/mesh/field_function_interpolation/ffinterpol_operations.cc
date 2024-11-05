// SPDX-FileCopyrightText: 2024 The OpenSn Authors <https://open-sn.github.io/opensn/>
// SPDX-License-Identifier: MIT

#include "lua/framework/lua.h"
#include "lua/framework/mesh/field_function_interpolation/ffinterpol.h"
#include "lua/framework/console/console.h"
#include "framework/field_functions/interpolation/ffinterpolation.h"
#include "framework/runtime.h"

namespace opensnlua
{

RegisterLuaFunctionInNamespace(FFInterpolationInitialize, fieldfunc, Initialize);
RegisterLuaFunctionInNamespace(FFInterpolationExecute, fieldfunc, Execute);

int
FFInterpolationInitialize(lua_State* L)
{
  const std::string fname = "fieldfunc.Initialize";
  LuaCheckArgs<size_t>(L, fname);

  // Get handle to field function
  const auto ffihandle = LuaArg<size_t>(L, 1);

  auto p_ffi = opensn::StackItemPtr(opensn::field_func_interpolation_stack, ffihandle, fname);

  p_ffi->Initialize();
  return LuaReturn(L);
}

int
FFInterpolationExecute(lua_State* L)
{
  const std::string fname = "fieldfunc.Execute";
  LuaCheckArgs<size_t>(L, fname);

  // Get handle to field function
  const auto ffihandle = LuaArg<size_t>(L, 1);

  auto p_ffi = opensn::StackItemPtr(opensn::field_func_interpolation_stack, ffihandle, fname);

  p_ffi->Execute();
  return LuaReturn(L);
}

} // namespace opensnlua
