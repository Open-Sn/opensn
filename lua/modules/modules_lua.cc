#include "modules_lua.h"
#include "framework/object_factory.h"
#include "lua/framework/console/console.h"
#include "cfem_diffusion/lua_cfem_diff_solver.h"
#include "dfem_diffusion/lua_dfem_diff_solver.h"
#include "fv_diffusion/lua_fv_diff_solver.h"
#include "mg_diffusion/lua_mg_diff_solver.h"
#include "linear_bolzmann_solvers/lbs_solver/lbs_lua_utils.h"
#include "config.h"

using namespace opensn;

namespace opensnlua
{

void
LoadRegisteredLuaItems()
{
  auto& L = console.GetConsoleState();

  luaL_openlibs(L);

  // Register version
  lua_pushstring(L, PROJECT_VERSION);
  lua_setglobal(L, "version");
  lua_pushinteger(L, PROJECT_MAJOR_VERSION);
  lua_setglobal(L, "major_version");
  lua_pushinteger(L, PROJECT_MINOR_VERSION);
  lua_setglobal(L, "minor_version");
  lua_pushinteger(L, PROJECT_PATCH_VERSION);
  lua_setglobal(L, "patch_version");

  // Registering functions
  RegisterLuaEntities(L);

  // Registering static-registration lua functions
  for (const auto& [key, entry] : console.GetLuaFunctionRegistry())
    Console::SetLuaFuncNamespaceTableStructure(key, entry.function_ptr);

  // Registering LuaFunctionWrappers
  for (const auto& [key, entry] : console.GetFunctionWrapperRegistry())
    if (entry.call_func)
      Console::SetLuaFuncWrapperNamespaceTableStructure(key);

  for (const auto& [key, value] : console.GetLuaConstantsRegistry())
    Console::SetLuaConstant(key, value);

  // Registering solver-function
  //                                    scope resolution tables
  const auto& object_maker = ObjectFactory::GetInstance();
  for (const auto& entry : object_maker.Registry())
    Console::SetObjectNamespaceTableStructure(entry.first);
}

void
RegisterLuaEntities(lua_State* L)
{
  opensnlua::lbs::RegisterLuaEntities(L);
}

} // namespace opensnlua
