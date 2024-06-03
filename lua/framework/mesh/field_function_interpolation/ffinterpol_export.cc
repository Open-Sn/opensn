// SPDX-FileCopyrightText: 2024 The OpenSn Authors <https://open-sn.github.io/opensn/>
// SPDX-License-Identifier: MIT

#include "framework/lua.h"
#include "lua/framework/mesh/field_function_interpolation/ffinterpol.h"
#include "framework/field_functions/interpolation/ffinterpolation.h"
#include "framework/console/console.h"
#include "framework/logging/log.h"
#include "framework/runtime.h"

namespace opensnlua
{

RegisterLuaFunctionInNamespace(FFInterpolationExportToCSV, fieldfunc, ExportToCSV);
RegisterLuaFunctionInNamespace(FFInterpolationExportToPython, fieldfunc, ExportToPython);

int
FFInterpolationExportToCSV(lua_State* L)
{
  // Get handle to field function
  const auto ffihandle = LuaArg<size_t>(L, 1);
  auto p_ffi = opensn::GetStackItemPtr(
    opensn::field_func_interpolation_stack, ffihandle, "fieldfunc.ExportToCSV");
  if (p_ffi->Type() != opensn::FieldFunctionInterpolationType::LINE)
    opensn::log.Log0Warning() << "ExportToCSV is only supported for LINE interpolators.";

  LuaCheckArgs<int>(L, "fieldfunc.ExportToCSV");
  auto base_name = LuaArgOptional<std::string>(L, 2, opensn::input_path.stem());
  p_ffi->ExportToCSV(base_name);

  return LuaReturn(L);
}

int
FFInterpolationExportToPython(lua_State* L)
{
  // Get handle to field function
  const auto ffihandle = LuaArg<size_t>(L, 1);
  auto p_ffi = opensn::GetStackItemPtr(
    opensn::field_func_interpolation_stack, ffihandle, "fieldfunc.ExportToPython");
  if (p_ffi->Type() != opensn::FieldFunctionInterpolationType::SLICE)
    opensn::log.Log0Warning() << "ExportToPython is only supported for SLICEinterpolators.";

  LuaCheckArgs<int>(L, "fieldfunc.ExportToPython");
  auto base_name = LuaArgOptional<std::string>(L, 2, opensn::input_path.stem());
  p_ffi->ExportToPython(base_name);

  return LuaReturn(L);
}

} // namespace opensnlua
