#include "framework/lua.h"
#include "unpartition_mesh_lua_utils.h"

#include "framework/runtime.h"
#include "framework/logging/log.h"
#include "framework/console/console.h"

namespace opensnlua
{

RegisterLuaFunctionAsIs(UnpartitionedMeshFinalizeEmpty);

int
UnpartitionedMeshFinalizeEmpty(lua_State* L)
{
  const std::string fname = __FUNCTION__;
  const int num_args = lua_gettop(L);
  if (num_args != 1)
    LuaPostArgAmountError(fname, 1, num_args);

  LuaCheckNilValue(fname, L, 1);

  const int handle = lua_tointeger(L, 1);

  auto& mesh =
    opensn::GetStackItem<opensn::UnpartitionedMesh>(opensn::unpartitionedmesh_stack, handle, fname);

  mesh.ComputeCentroidsAndCheckQuality();
  mesh.BuildMeshConnectivity();

  return 0;
}

} // namespace opensnlua
