// SPDX-FileCopyrightText: 2024 The OpenSn Authors <https://open-sn.github.io/opensn/>
// SPDX-License-Identifier: MIT

#include "framework/lua.h"
#include "framework/field_functions/interpolation/ffinter_point.h"
#include "framework/field_functions/interpolation/ffinter_line.h"
#include "framework/field_functions/interpolation/ffinter_volume.h"
#include "framework/runtime.h"
#include "framework/logging/log.h"
#include "ffinterpol_lua.h"
#include "framework/console/console.h"

using namespace opensn;

namespace opensnlua
{

RegisterLuaFunctionInNamespace(FFInterpolationGetValue, fieldfunc, GetValue);

int
FFInterpolationGetValue(lua_State* L)
{
  const std::string fname = "fieldfunc.GetValue";
  LuaCheckArgs<size_t>(L, fname);

  // Get handle to field function
  const auto ffihandle = LuaArg<size_t>(L, 1);

  auto p_ffi = opensn::GetStackItemPtr(opensn::field_func_interpolation_stack, ffihandle, fname);

  if (p_ffi->Type() == FieldFunctionInterpolationType::POINT)
  {
    auto& cur_ffi_point = dynamic_cast<opensn::FieldFunctionInterpolationPoint&>(*p_ffi);
    double value = cur_ffi_point.GetPointValue();
    return LuaReturn(L, value);
  }
  else if (p_ffi->Type() == FieldFunctionInterpolationType::LINE)
  {
    auto& cur_ffi_line = dynamic_cast<opensn::FieldFunctionInterpolationLine&>(*p_ffi);

    std::vector<std::vector<double>> vals;
    vals.resize(cur_ffi_line.GetFieldFunctions().size());
    for (int ff = 0; ff < cur_ffi_line.GetFieldFunctions().size(); ff++)
    {
      const auto& ff_ctx = cur_ffi_line.GetFFContexts()[ff];
      vals[ff] = ff_ctx.interpolation_points_values;
    }
    return LuaReturn(L, vals);
  }
  else if (p_ffi->Type() == FieldFunctionInterpolationType::VOLUME)
  {
    auto& cur_ffi_volume = dynamic_cast<opensn::FieldFunctionInterpolationVolume&>(*p_ffi);
    double value = cur_ffi_volume.GetOpValue();
    return LuaReturn(L, value);
  }
  else
  {
    opensn::log.Log0Warning()
      << fname + " is currently only supported for POINT, LINE and VOLUME interpolator types.";
  }

  return LuaReturn(L);
}

} // namespace opensnlua
