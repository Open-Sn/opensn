// SPDX-FileCopyrightText: 2024 The OpenSn Authors <https://open-sn.github.io/opensn/>
// SPDX-License-Identifier: MIT

#include "framework/lua.h"
#include "framework/field_functions/field_function_grid_based.h"
#include "framework/runtime.h"
#include "framework/logging/log.h"
#include "field_functions_lua.h"
#include "framework/console/console.h"

using namespace opensn;

namespace opensnlua
{

RegisterLuaFunctionNamespace(ExportFieldFunctionToVTK, fieldfunc, ExportToVTK);
RegisterLuaFunctionNamespace(ExportMultiFieldFunctionToVTK, fieldfunc, ExportToVTKMulti);

int
ExportFieldFunctionToVTK(lua_State* L)
{
  const std::string fname = "fieldfunc.ExportToVTK";
  LuaCheckArgs<size_t, std::string>(L, fname);

  auto ff_handle = LuaArg<size_t>(L, 1);
  auto base_name = LuaArg<std::string>(L, 2);

  auto ff_base = opensn::GetStackItemPtr(opensn::field_function_stack, ff_handle, fname);
  auto ff = std::dynamic_pointer_cast<FieldFunctionGridBased>(ff_base);

  OpenSnLogicalErrorIf(not ff, "Only grid-based field functions can be exported");

  //  ff->ExportToVTK(base_name);
  FieldFunctionGridBased::ExportMultipleToVTK(base_name, {ff});

  return LuaReturn(L);
}

int
ExportMultiFieldFunctionToVTK(lua_State* L)
{
  const std::string fname = "fieldfunc.ExportToVTKMulti";
  LuaCheckArgs<std::vector<size_t>, std::string>(L, fname);

  auto ff_handles = LuaArg<std::vector<size_t>>(L, 1);
  auto base_name = LuaArg<std::string>(L, 2);

  std::vector<std::shared_ptr<const FieldFunctionGridBased>> ffs;
  ffs.reserve(ff_handles.size());
  for (std::size_t i = 0; i < ff_handles.size(); ++i)
  {
    std::shared_ptr<FieldFunction> ff_base =
      opensn::GetStackItemPtr(opensn::field_function_stack, ff_handles[i], fname);
    auto ff = std::dynamic_pointer_cast<FieldFunctionGridBased>(ff_base);
    OpenSnLogicalErrorIf(not ff, "Only grid-based field functions can be exported");

    ffs.push_back(ff);
  }

  FieldFunctionGridBased::ExportMultipleToVTK(base_name, ffs);

  return LuaReturn(L);
}

} // namespace opensnlua
