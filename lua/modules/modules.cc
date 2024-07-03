// SPDX-FileCopyrightText: 2024 The OpenSn Authors <https://open-sn.github.io/opensn/>
// SPDX-License-Identifier: MIT

#include "lua/modules/modules.h"
#include "framework/object_factory.h"
#include "lua/framework/lua.h"
#include "lua/framework/console/console.h"
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
  LuaSetGlobal(L, "version", PROJECT_VERSION);
  LuaSetGlobal(L, "major_version", PROJECT_MAJOR_VERSION);
  LuaSetGlobal(L, "minor_version", PROJECT_MINOR_VERSION);
  LuaSetGlobal(L, "patch_version", PROJECT_PATCH_VERSION);

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
  RegisterLuaEntitiesLBS(L);
}

} // namespace opensnlua
