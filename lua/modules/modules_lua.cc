#include "modules_lua.h"
#include "framework/object_factory.h"

#include "diffusion_solver/diffusion_lua.h"
#include "cfem_diffusion/ds_lua_utils.h"
#include "dfem_diffusion/ip_lua_utils.h"
#include "fv_diffusion/ds_lua_utils.h"
#include "mg_diffusion/mgds_lua_utils.h"
#include "linear_bolzmann_solvers/lbs_solver/lbs_lua_utils.h"
#include "config.h"

using namespace opensn;

namespace opensnlua
{

void
LoadRegisteredLuaItems()
{
  // Initializing console
  auto console = opensn::Chi::console;

  auto& L = console.GetConsoleState();

  luaL_openlibs(L);

  // Register version
  lua_pushstring(L, PROJECT_VERSION);
  lua_setglobal(L, "chi_version");
  lua_pushinteger(L, PROJECT_MAJOR_VERSION);
  lua_setglobal(L, "chi_major_version");
  lua_pushinteger(L, PROJECT_MINOR_VERSION);
  lua_setglobal(L, "chi_minor_version");
  lua_pushinteger(L, PROJECT_PATCH_VERSION);
  lua_setglobal(L, "chi_patch_version");

  // Registering functions
  RegisterLuaEntities(L);

  // Registering static-registration lua functions
  for (const auto& [key, entry] : console.GetLuaFunctionRegistry())
    Console::SetLuaFuncNamespaceTableStructure(key, entry.function_ptr);

  // Registering LuaFunctionWrappers
  for (const auto& [key, entry] : console.GetFunctionWrapperRegistry())
    if (entry.call_func) Console::SetLuaFuncWrapperNamespaceTableStructure(key);

  for (const auto& [key, value] : console.GetLuaConstantsRegistry())
    Console::SetLuaConstant(key, value);

  // Registering solver-function
  //                                    scope resolution tables
  const auto& object_maker = ChiObjectFactory::GetInstance();
  for (const auto& entry : object_maker.Registry())
    Console::SetObjectNamespaceTableStructure(entry.first);
}

void
RegisterLuaEntities(lua_State* L)
{
  opensnlua::lbs::RegisterLuaEntities(L);
  opensnlua::diffusion_solver::RegisterLuaEntities(L);

  opensnlua::cfem_diffusion::RegisterLuaEntities(L);
  opensnlua::dfem_diffusion::RegisterLuaEntities(L);

  opensnlua::mg_diffusion::RegisterLuaEntities(L);
  opensnlua::fv_diffusion::RegisterLuaEntities(L);
}

} // namespace opensnlua
