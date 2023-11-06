#include "framework/chi_lua.h"

#include "framework/mesh/MeshGenerator/MeshGenerator.h"

#include "framework/console/chi_console.h"

#include "framework/chi_runtime.h"

namespace chi_mesh::lua_utils
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

  auto& generator = Chi::GetStackItem<MeshGenerator>(Chi::object_stack, handle, fname);
  generator.Execute();

  return 0;
}

} // namespace chi_mesh::lua_utils
