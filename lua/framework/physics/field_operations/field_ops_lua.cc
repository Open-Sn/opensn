#include "field_ops_lua.h"
#include "framework/physics/field_operations/field_operation.h"

#include "framework/console/console.h"

using namespace opensn;

namespace opensnlua
{

RegisterLuaFunctionAsIs(chiFieldOperationExecute);

int
chiFieldOperationExecute(lua_State* L)
{
  const std::string fname = __FUNCTION__;
  const int num_args = lua_gettop(L);
  if (num_args != 1) LuaPostArgAmountError(fname, 1, num_args);

  LuaCheckNilValue(fname, L, 1);

  const size_t handle = lua_tointeger(L, 1);

  auto& operation =
    opensn::Chi::GetStackItem<FieldOperation>(opensn::Chi::object_stack, handle, fname);

  operation.Execute();

  return 0;
}

} // namespace opensnlua
