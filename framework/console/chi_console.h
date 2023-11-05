#pragma once

extern "C"
{
#include <lua.h>
#include <lualib.h>
#include <lauxlib.h>
}

#include "opensn/framework/parameters/parameter_block.h"
#include "opensn/framework/parameters/input_parameters.h"
#include "opensn/framework/logging/chi_log_exceptions.h"

#include <vector>
#include <string>
#include <map>
#include <stack>

class Chi;

#ifdef OPENSN_WITH_LUA

/**Small utility macro for joining two words.*/
#define ConsoleJoinWordsA(x, y) x##y
/**IDK why this is needed. Seems like counter doesnt work properly without it*/
#define ConsoleJoinWordsB(x, y) ConsoleJoinWordsA(x, y)

/**Macro for registering a lua_CFunction within the Console
 * singleton, with the function being in the global namespace. Example:
 * \code
 * ConsoleRegisterLuaFunction(chiSolverInitialize);
 * \endcode
 *
 * \note Remember to include the header "console/chi_console.h".
 * The name supplied to this function cannot have scope resolution operators,
 * i.e., "::".*/
#define RegisterLuaFunctionAsIs(func_name)                                                         \
  static char ConsoleJoinWordsB(unique_var_name_luacfunc_##func_name##_, __COUNTER__) =            \
    chi::Console::AddFunctionToRegistryGlobalNamespace(#func_name, func_name)

/**Macro for registering a lua_CFunction within the Console
* singleton.
\param function LuaCFunction. The function to use.
\param namespace_name NonQuotedString. May include scope resolution
\param func_name NonQuotedString. The name of the function as it will appear in
                 the lua console.
*/
#define RegisterLuaFunction(function, namespace_name, func_name)                                   \
  static char ConsoleJoinWordsB(unique_var_name_luacfunc_##func_name##_, __COUNTER__) =            \
    chi::Console::AddFunctionToRegistryInNamespaceWithName(function, #namespace_name, #func_name)

#define RegisterWrapperFunction(namespace_name, name_in_lua, syntax_function, actual_function)     \
  static char ConsoleJoinWordsB(unique_var_name_luacfunc_##name_in_lua##_, __COUNTER__) =          \
    chi::Console::AddWrapperToRegistryInNamespaceWithName(                                         \
      #namespace_name, #name_in_lua, syntax_function, actual_function)

#define RegisterLuaConstant(namespace_name, name_in_lua, value)                                    \
  static char ConsoleJoinWordsB(unique_var_name_luaconst_##namespace_name##_##name_in_lua,         \
                                __COUNTER__) =                                                     \
    chi::Console::AddLuaConstantToRegistry(#namespace_name, #name_in_lua, value)

#define RegisterLuaConstantAsIs(name_in_lua, value)                                                \
  static char ConsoleJoinWordsB(unique_var_name_luaconst_##name_in_lua, __COUNTER__) =             \
    chi::Console::AddLuaConstantToRegistry("", #name_in_lua, value)

#else
#define RegisterLuaFunctionAsIs(func_name)
#define RegisterLuaFunction(function, namespace_name, func_name)
#define RegisterWrapperFunction(namespace_name, name_in_lua, syntax_function, actual_function)
#define RegisterLuaConstant(namespace_name, name_in_lua, value)
#define RegisterLuaConstantAsIs(name_in_lua, value)
#endif

namespace chi_physics
{
class Solver;
}
namespace chi
{
struct RegistryStatuses;
}

namespace chi
{

/**
 * Class for handling the console and scripting.
 */
class Console
{
public:
  using WrapperGetInParamsFunc = chi::InputParameters (*)();
  using WrapperCallFunc = chi::ParameterBlock (*)(const chi::InputParameters&);

private:
#ifdef OPENSN_WITH_LUA
  struct LuaFunctionRegistryEntry
  {
    lua_CFunction function_ptr;
    std::string function_raw_name;
  };
  struct LuaFuncWrapperRegEntry
  {
    WrapperGetInParamsFunc get_in_params_func = nullptr;
    WrapperCallFunc call_func = nullptr;
  };

  /// Pointer to lua console state
  lua_State* console_state_;
#endif
  /// Buffer of commands to execute
  std::vector<std::string> command_buffer_;
  static Console instance_;

#ifdef OPENSN_WITH_LUA
  std::map<std::string, LuaFunctionRegistryEntry> lua_function_registry_;

  std::map<std::string, LuaFuncWrapperRegEntry> function_wrapper_registry_;
#endif

  std::map<std::string, chi_data_types::Varying> lua_constants_registry_;

  Console() noexcept;

private:
  friend class ::Chi;

public:
  /**
   * Access to the singleton
   */
  static Console& GetInstance() noexcept;

#ifdef OPENSN_WITH_LUA
  lua_State*& GetConsoleState() { return console_state_; }
#endif

  std::vector<std::string>& GetCommandBuffer() { return command_buffer_; }

#ifdef OPENSN_WITH_LUA
  const std::map<std::string, LuaFunctionRegistryEntry>& GetLuaFunctionRegistry() const
  {
    return lua_function_registry_;
  }

  const std::map<std::string, LuaFuncWrapperRegEntry>& GetFunctionWrapperRegistry() const
  {
    return function_wrapper_registry_;
  }

  const std::map<std::string, chi_data_types::Varying> & GetLuaConstantsRegistry() const
  {
    return lua_constants_registry_;
  }
#endif

  /**
   * Executes the loop for the console.
   */
  void RunConsoleLoop(char* fileName = nullptr) const;
  int ExecuteFile(const std::string& fileName, int argc, char** argv) const;
  void PostMPIInfo(int location_id, int number_of_processes) const;

private:
  /**
   * Basic addition to registry. Used by the other public methods to registry a text-key to a lua
   * function.
   */
  static void AddFunctionToRegistry(const std::string& name_in_lua, lua_CFunction function_ptr);

public:
#ifdef OPENSN_WITH_LUA
  /**
   * Adds a lua_CFunction to the registry.
   */
  static char AddFunctionToRegistryGlobalNamespace(const std::string& raw_name_in_lua,
                                                   lua_CFunction function_ptr);

  /**
   * Adds a lua_CFunction to the registry. With namespace-table analogy.
   */
  static char AddFunctionToRegistryInNamespaceWithName(lua_CFunction function_ptr,
                                                       const std::string& namespace_name,
                                                       const std::string& function_name);

  /**
   * Adds a constant to the lua state.
   */
  static char AddLuaConstantToRegistry(const std::string& namespace_name,
                                       const std::string& constant_name,
                                       const chi_data_types::Varying& value);
#endif

  /**
   * A default function for returning empty input parameters.
   */
  static InputParameters DefaultGetInParamsFunc();

#ifdef OPENSN_WITH_LUA
  /**
   * Adds a function wrapper to the lua registry.
   */
  static char AddWrapperToRegistryInNamespaceWithName(const std::string& namespace_name,
                                                      const std::string& name_in_lua,
                                                      WrapperGetInParamsFunc syntax_function,
                                                      WrapperCallFunc actual_function,
                                                      bool ignore_null_call_func = false);

  /**
   * Formats a namespace structure as table.
   */
  static void SetLuaFuncNamespaceTableStructure(const std::string& full_lua_name,
                                                lua_CFunction function_ptr);

  /**
   * Formats a namespace structure as a table, but the last entry is a function call.
   */
  static void SetLuaFuncWrapperNamespaceTableStructure(const std::string& full_lua_name);

  /**
   * Formats a namespace structure as a table, but the last entry contains a "Create" function and
   * a type.
   */
  static void SetObjectNamespaceTableStructure(const std::string& full_lua_name);

  /**
   * Makes sure a table structure exists for the list of table names.
   */
  static void FleshOutLuaTableStructure(const std::vector<std::string>& table_names);

  /**
   * Sets a lua constant in the lua state.
   */
  static void SetLuaConstant(const std::string& constant_name,
                             const chi_data_types::Varying& value);
#endif

  /**
   * Flushes any commands in the command buffer.
   */
  void FlushConsole();

#ifdef OPENSN_WITH_LUA
  /**
   * Generic entry point for wrapper calls.
   */
  static int LuaWrapperCall(lua_State* L);

  /**
   * Dumps the object registry to stdout.
   */
  void DumpRegister() const;

  /**
   * Given an old status, will update the bindings for only newly registered items.
   */
  void UpdateConsoleBindings(const chi::RegistryStatuses& old_statuses);
#endif
};

} // namespace chi
