// SPDX-FileCopyrightText: 2024 The OpenSn Authors <https://open-sn.github.io/opensn/>
// SPDX-License-Identifier: MIT

#include "lua/framework/mesh/logical_volume/logic_volume.h"
#include "framework/mesh/logical_volume/logical_volume.h"
#include "framework/runtime.h"
#include "lua/framework/console/console.h"

using namespace opensn;

namespace opensnlua
{

RegisterLuaFunctionInNamespace(LogVolPointSense, logvol, PointSense);

int
LogVolPointSense(lua_State* L)
{
  const std::string fname = "logvol.PointSense";
  LuaCheckArgs<int, Vector3>(L, fname);

  const auto lv_handle = LuaArg<int>(L, 1);
  auto point = LuaArg<Vector3>(L, 2);

  const auto& lv = opensn::GetStackItem<LogicalVolume>(opensn::object_stack, lv_handle, fname);

  auto ret_val = lv.Inside(point);
  return LuaReturn(L, ret_val);
}

} // namespace opensnlua
