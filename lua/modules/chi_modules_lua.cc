#include "chi_modules_lua.h"
#include "framework/ChiObjectFactory.h"

#include "diffusion_solver/diffusion_lua.h"
#include "cfem_diffusion/ds_lua_utils.h"
#include "dfem_diffusion/ip_lua_utils.h"
#include "fv_diffusion/ds_lua_utils.h"
#include "mg_diffusion/mgds_lua_utils.h"
#include "linear_bolzmann_solvers/lbs_solver/lbs_lua_utils.h"
#include "config.h"

void
chi_modules::lua_utils::LoadRegisteredLuaItems()
{
  // Initializing console
  auto console = Chi::console;

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
  chi_modules::lua_utils::RegisterLuaEntities(L);

  // Registering static-registration lua functions
  for (const auto& [key, entry] : console.GetLuaFunctionRegistry())
    chi::Console::SetLuaFuncNamespaceTableStructure(key, entry.function_ptr);

  // Registering LuaFunctionWrappers
  for (const auto& [key, entry] : console.GetFunctionWrapperRegistry())
    if (entry.call_func) chi::Console::SetLuaFuncWrapperNamespaceTableStructure(key);

  for (const auto& [key, value] : console.GetLuaConstantsRegistry())
    chi::Console::SetLuaConstant(key, value);

  // Registering solver-function
  //                                    scope resolution tables
  const auto& object_maker = ChiObjectFactory::GetInstance();
  for (const auto& entry : object_maker.Registry())
    chi::Console::SetObjectNamespaceTableStructure(entry.first);
}

void
chi_modules::lua_utils::RegisterLuaEntities(lua_State* L)
{
  lbs::common_lua_utils::RegisterLuaEntities(L);
  diffusion_solver::lua_utils::RegisterLuaEntities(L);

  cfem_diffusion::cfem_diffusion_lua_utils::RegisterLuaEntities(L);
  dfem_diffusion::dfem_diffusion_lua_utils::RegisterLuaEntities(L);

  mg_diffusion::mgd_lua_utils::RegisterLuaEntities(L);
  fv_diffusion::fv_diffusion_lua_utils::RegisterLuaEntities(L);
}
