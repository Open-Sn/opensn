#include "framework/lua.h"

#include "framework/mesh/mesh_generator/mesh_generator.h"

#include "framework/console/console.h"

#include "framework/runtime.h"

using namespace opensn;

namespace opensnlua
{

int chiMeshGeneratorExecute(lua_State* L);

RegisterLuaFunction(chiMeshGeneratorExecute, chi_mesh::MeshGenerator, Execute);

int
chiMeshGeneratorExecute(lua_State* L)
{
  const std::string fname = __FUNCTION__;
  const int num_args = lua_gettop(L);
  if (num_args != 1) LuaPostArgAmountError(fname, 1, num_args);

  LuaCheckNilValue(fname, L, 1);
  LuaCheckIntegerValue(fname, L, 1);

  const size_t handle = lua_tointeger(L, 1);

  auto& generator = opensn::GetStackItem<MeshGenerator>(opensn::object_stack, handle, fname);
  generator.Execute();

  return 0;
}

} // namespace opensnlua
