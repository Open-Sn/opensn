#pragma once

extern "C"
{
#include <lua.h>
#include <lualib.h>
#include <lauxlib.h>
}
#include "chi_console_structs.h"
#include "parameters/parameter_block.h"
#include "parameters/input_parameters.h"

#include "chi_log_exceptions.h"

#include <vector>
#include <string>
#include <map>
#include <stack>

class Chi;

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

namespace chi_physics
{
class Solver;
}
namespace chi
{
struct RegistryStatuses;
}

// #############################################################################
// CLASS DEF
namespace chi
{

/** Class for handling the console and scripting.*/
class Console
{
public:
  using WrapperGetInParamsFunc = chi::InputParameters (*)();
  using WrapperCallFunc = chi::ParameterBlock (*)(const chi::InputParameters&);

private:
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

private:
  lua_State* console_state_; ///< Pointer to lua console state

  std::vector<std::string> command_buffer_; ///< Buffer of commands to execute
  static Console instance_;

  std::map<std::string, LuaFunctionRegistryEntry> lua_function_registry_;

  std::map<std::string, LuaFuncWrapperRegEntry> function_wrapper_registry_;

  std::map<std::string, chi_data_types::Varying> lua_constants_registry_;

  // 00
  Console() noexcept;

private:
  friend class ::Chi;
  void LoadRegisteredLuaItems();

public:
  static Console& GetInstance() noexcept;

  lua_State*& GetConsoleState() { return console_state_; }
  std::vector<std::string>& GetCommandBuffer() { return command_buffer_; }

  const std::map<std::string, LuaFunctionRegistryEntry>& GetLuaFunctionRegistry() const
  {
    return lua_function_registry_;
  }

  const std::map<std::string, LuaFunctionRegistryEntry>& GetFunctionWrapperRegistry() const
  {
    return lua_function_registry_;
  }

  // 01 Loop
  void RunConsoleLoop(char* fileName = nullptr) const;
  // 02 Utilities
  int ExecuteFile(const std::string& fileName, int argc, char** argv) const;
  void PostMPIInfo(int location_id, int number_of_processes) const;

private:
  /**Basic addition to registry. Used by the other public methods
   * to registry a text-key to a lua function.*/
  static void AddFunctionToRegistry(const std::string& name_in_lua, lua_CFunction function_ptr);

public:
  /**\brief Adds a lua_CFunction to the registry.*/
  static char AddFunctionToRegistryGlobalNamespace(const std::string& raw_name_in_lua,
                                                   lua_CFunction function_ptr);

  /**\brief Adds a lua_CFunction to the registry. With namespace-table
   * analogy.*/
  static char AddFunctionToRegistryInNamespaceWithName(lua_CFunction function_ptr,
                                                       const std::string& namespace_name,
                                                       const std::string& function_name);

  /**\brief Adds a constant to the lua state.*/
  static char AddLuaConstantToRegistry(const std::string& namespace_name,
                                       const std::string& constant_name,
                                       const chi_data_types::Varying& value);

  /**\brief A default function for returning empty input parameters. */
  static InputParameters DefaultGetInParamsFunc();

  /**\brief Adds a function wrapper to the lua registry.*/
  static char AddWrapperToRegistryInNamespaceWithName(const std::string& namespace_name,
                                                      const std::string& name_in_lua,
                                                      WrapperGetInParamsFunc syntax_function,
                                                      WrapperCallFunc actual_function,
                                                      bool ignore_null_call_func = false);

  /**\brief Formats a namespace structure as table.*/
  static void SetLuaFuncNamespaceTableStructure(const std::string& full_lua_name,
                                                lua_CFunction function_ptr);

  /**\brief Formats a namespace structure as a table, but the last entry
   * is a function call.*/
  static void SetLuaFuncWrapperNamespaceTableStructure(const std::string& full_lua_name);

  /**\brief Formats a namespace structure as a table, but the last entry
   * contains a "Create" function and a type.*/
  static void SetObjectNamespaceTableStructure(const std::string& full_lua_name);

  /**\brief Makes sure a table structure exists for the list of table names.*/
  static void FleshOutLuaTableStructure(const std::vector<std::string>& table_names);

  /**Sets a lua constant in the lua state.*/
  static void SetLuaConstant(const std::string& constant_name,
                             const chi_data_types::Varying& value);

  // 03
  /**\brief Flushes any commands in the command buffer.*/
  void FlushConsole();
  // 05 Memory
  /**\brief Get current memory usage.*/
  static CSTMemory GetMemoryUsage();
  /**\brief Get current memory usage in Megabytes.*/
  static double GetMemoryUsageInMB();

  // fwrapper_call
  /**\brief Generic entry point for wrapper calls.*/
  static int LuaWrapperCall(lua_State* L);

  /**\brief Dumps the object registry to stdout.*/
  void DumpRegister() const;

  /**Given an old status, will update the bindings for only newly registered
   * items.*/
  void UpdateConsoleBindings(const chi::RegistryStatuses& old_statuses);
};
} // namespace chi


