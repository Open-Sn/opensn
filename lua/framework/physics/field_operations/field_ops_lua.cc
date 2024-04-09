// SPDX-FileCopyrightText: 2024 The OpenSn Authors <https://open-sn.github.io/opensn/>
// SPDX-License-Identifier: MIT

#include "field_ops_lua.h"
#include "framework/field_functions/operations/field_operation.h"
#include "framework/console/console.h"

using namespace opensn;

namespace opensnlua
{

RegisterLuaFunctionNamespace(FieldOperationExecute, fieldfunc, FieldOperationExecute);

int
FieldOperationExecute(lua_State* L)
{
  const std::string fname = "fieldfunc.FieldOperationExecute";
  LuaCheckArgs<size_t>(L, fname);

  const auto handle = LuaArg<size_t>(L, 1);

  auto& operation = opensn::GetStackItem<FieldOperation>(opensn::object_stack, handle, fname);

  operation.Execute();

  return LuaReturn(L);
}

} // namespace opensnlua
