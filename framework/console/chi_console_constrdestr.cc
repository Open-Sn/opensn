#include "console/chi_console.h"

#include "lua/chi_modules_lua.h"

#include "chi_configuration.h"

#include "ChiObjectFactory.h"

//###################################################################
/**Access to the singleton*/
chi::Console& chi::Console::GetInstance() noexcept
{
  static Console singleton;
  return singleton;
}

//###################################################################
/** Default constructor for the console*/
chi::Console::Console() noexcept :
  console_state_(luaL_newstate())
{

}

//###################################################################
/**Registers all lua items so that they are available in the console.*/
void chi::Console::LoadRegisteredLuaItems()
{
  //=================================== Initializing console
  auto& L = GetConsoleState();

  luaL_openlibs(L);

  //=================================== Register version
  lua_pushstring(L, PROJECT_VERSION);      lua_setglobal(L,"chi_version");
  lua_pushinteger(L,PROJECT_MAJOR_VERSION);lua_setglobal(L,"chi_major_version");
  lua_pushinteger(L,PROJECT_MINOR_VERSION);lua_setglobal(L,"chi_minor_version");
  lua_pushinteger(L,PROJECT_PATCH_VERSION);lua_setglobal(L,"chi_patch_version");

  //=================================== Registering functions
  chi_modules::lua_utils::RegisterLuaEntities(L);

  //=================================== Registering static-registration
  //                                    lua functions
  for (const auto& [key, entry] : lua_function_registry_)
    SetLuaFuncNamespaceTableStructure(key, entry.function_ptr);

  //=================================== Registering LuaFunctionWrappers
  for (const auto& [key, entry] : function_wrapper_registry_)
    if (entry.call_func)
      SetLuaFuncWrapperNamespaceTableStructure(key);

  for (const auto& [key, value] : lua_constants_registry_)
    SetLuaConstant(key, value);

  //=================================== Registering solver-function
  //                                    scope resolution tables
  const auto& object_maker = ChiObjectFactory::GetInstance();
  for (const auto& entry : object_maker.Registry())
    SetObjectNamespaceTableStructure(entry.first);

}


