#pragma once

extern "C"
{
#include <lua.h>
#include <lualib.h>
#include <lauxlib.h>
}

#include "framework/parameters/parameter_block.h"
#include "framework/parameters/input_parameters.h"
#include "framework/logging/log_exceptions.h"

#include <vector>
#include <string>
#include <map>
#include <stack>

/**Small utility macro for joining two words.*/
#define ConsoleJoinWordsA(x, y) x##y
/**IDK why this is needed. Seems like counter doesnt work properly without it*/
#define ConsoleJoinWordsB(x, y) ConsoleJoinWordsA(x, y)

/**Macro for registering a lua_CFunction within the Console
 * singleton, with the function being in the global namespace. Example:
 * \code
 * ConsoleRegisterLuaFunction(SolverInitialize);
 * \endcode
 *
 * \note Remember to include the header "console/chi_console.h".
 * The name supplied to this function cannot have scope resolution operators,
 * i.e., "::".*/
#define RegisterLuaFunctionAsIs(func_name)                                                         \
  static char ConsoleJoinWordsB(unique_var_name_luacfunc_##func_name##_, __COUNTER__) =            \
    opensnlua::Console::AddFunctionToRegistryGlobalNamespace(#func_name, func_name)

/**Macro for registering a lua_CFunction within the Console
* singleton.
\param function LuaCFunction. The function to use.
\param namespace_name NonQuotedString. May include scope resolution
\param func_name NonQuotedString. The name of the function as it will appear in
                 the lua console.
*/
#define RegisterLuaFunctionNamespace(function, namespace_name, func_name)                          \
  static char ConsoleJoinWordsB(unique_var_name_luacfunc_##func_name##_, __COUNTER__) =            \
    opensnlua::Console::AddFunctionToRegistryInNamespaceWithName(                                  \
      function, #namespace_name, #func_name)

#define RegisterLuaFunction(function, func_name)                                                   \
  static char ConsoleJoinWordsB(unique_var_name_luacfunc_##func_name##_, __COUNTER__) =            \
    opensnlua::Console::AddFunctionToRegistryInNamespaceWithName(function, #func_name)

#define RegisterWrapperFunctionNamespace(                                                          \
  namespace_name, name_in_lua, syntax_function, actual_function)                                   \
  static char ConsoleJoinWordsB(unique_var_name_luacfunc_##name_in_lua##_, __COUNTER__) =          \
    opensnlua::Console::AddWrapperToRegistryInNamespaceWithName(                                   \
      #namespace_name, #name_in_lua, syntax_function, actual_function)

#define RegisterWrapperFunction(name_in_lua, syntax_function, actual_function)                     \
  static char ConsoleJoinWordsB(unique_var_name_luacfunc_##name_in_lua##_, __COUNTER__) =          \
    opensnlua::Console::AddWrapperToRegistryInNamespaceWithName(                                   \
      #name_in_lua, syntax_function, actual_function)

#define RegisterLuaConstant(namespace_name, name_in_lua, value)                                    \
  static char ConsoleJoinWordsB(unique_var_name_luaconst_##namespace_name##_##name_in_lua,         \
                                __COUNTER__) =                                                     \
    opensnlua::Console::AddLuaConstantToRegistry(#namespace_name, #name_in_lua, value)

#define RegisterLuaConstantAsIs(name_in_lua, value)                                                \
  static char ConsoleJoinWordsB(unique_var_name_luaconst_##name_in_lua, __COUNTER__) =             \
    opensnlua::Console::AddLuaConstantToRegistry("", #name_in_lua, value)

namespace opensn
{
class Solver;
}

namespace opensnlua
{

struct RegistryStatuses;

/**
 * Class for handling the console and scripting.
 */
class Console
{
public:
  using WrapperGetInParamsFunc = opensn::InputParameters (*)();
  using WrapperCallFunc = opensn::ParameterBlock (*)(const opensn::InputParameters&);

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

  /// Pointer to lua console state
  lua_State* console_state_;
  /// Buffer of commands to execute
  std::vector<std::string> command_buffer_;
  static Console instance_;

  std::map<std::string, LuaFunctionRegistryEntry> lua_function_registry_;

  std::map<std::string, LuaFuncWrapperRegEntry> function_wrapper_registry_;

  std::map<std::string, opensn::Varying> lua_constants_registry_;

  Console() noexcept;

public:
  /**
   * Access to the singleton
   */
  static Console& GetInstance() noexcept;

  lua_State*& GetConsoleState() { return console_state_; }

  std::vector<std::string>& GetCommandBuffer() { return command_buffer_; }

  const std::map<std::string, LuaFunctionRegistryEntry>& GetLuaFunctionRegistry() const
  {
    return lua_function_registry_;
  }

  const std::map<std::string, LuaFuncWrapperRegEntry>& GetFunctionWrapperRegistry() const
  {
    return function_wrapper_registry_;
  }

  const std::map<std::string, opensn::Varying>& GetLuaConstantsRegistry() const
  {
    return lua_constants_registry_;
  }

  /**
   * Executes the loop for the console.
   */
  void RunConsoleLoop(char* fileName = nullptr) const;

  /**
   * Executes the given file in the Lua engine.
   * \author Jan
   */
  int ExecuteFile(const std::string& fileName, int argc, char** argv) const;

  /**
   * Pushes location id and number of processes to lua state.
   */
  void PostMPIInfo(int location_id, int number_of_processes) const;

private:
  /**
   * Basic addition to registry. Used by the other public methods to registry a text-key to a lua
   * function.
   */
  static void AddFunctionToRegistry(const std::string& name_in_lua, lua_CFunction function_ptr);

public:
  /**
   * Adds a lua_CFunction to the registry. The registry of functions gets parsed into the lua
   * console when `chi::Initialize` is called. This particular function will strip the namespace
   * from the the parameter `raw_name_in_lua` and cause the function to be registered in the
   * global namespace of the lua console.
   */
  static char AddFunctionToRegistryGlobalNamespace(const std::string& raw_name_in_lua,
                                                   lua_CFunction function_ptr);

  /**
   * Adds a lua_CFunction to the registry. The registry of functions gets parsed into the lua
   * console when `chi::Initialize` is called. The full path of the function will be derived from
   * `namespace_name` + "::" + `function_name`.
   */
  static char AddFunctionToRegistryInNamespaceWithName(lua_CFunction function_ptr,
                                                       const std::string& namespace_name,
                                                       const std::string& function_name);

  /**
   * \brief Adds a constant to the lua state. Prepending the constant within a namespace is
   * optional.
   */
  static char AddLuaConstantToRegistry(const std::string& namespace_name,
                                       const std::string& constant_name,
                                       const opensn::Varying& value);

  /**
   * A default function for returning empty input parameters.
   */
  static opensn::InputParameters DefaultGetInParamsFunc();

  /**
   * Adds a function wrapper to the lua registry.
   */
  static char AddWrapperToRegistryInNamespaceWithName(const std::string& namespace_name,
                                                      const std::string& name_in_lua,
                                                      WrapperGetInParamsFunc syntax_function,
                                                      WrapperCallFunc actual_function,
                                                      bool ignore_null_call_func = false);

  static char AddWrapperToRegistryInNamespaceWithName(const std::string& name_in_lua,
                                                      WrapperGetInParamsFunc syntax_function,
                                                      WrapperCallFunc actual_function,
                                                      bool ignore_null_call_func = false);

  /**
   * Formats a namespace structure as table.
   */
  static void SetLuaFuncNamespaceTableStructure(const std::string& full_lua_name,
                                                lua_CFunction function_ptr);

  /**
   * Sets/Forms a table structure that mimics the namespace structure of a string. For example
   * the string "sing::sob::nook::Tigger" will be assigned a table structure
   * `sing.sob.nook.Tigger = "sing::sob::nook::Tigger"`. Then finally assigns lua call to this
   * table.
   */
  static void SetLuaFuncWrapperNamespaceTableStructure(const std::string& full_lua_name);

  /**
   * Sets/Forms a table structure that mimics the namespace structure of a string. For example
   * the string "sing::sob::nook::Tigger" will be assigned a table structure
   * `sing.sob.nook.Tigger = "sing::sob::nook::Tigger"`.
   */
  static void SetObjectNamespaceTableStructure(const std::string& full_lua_name);

  /**
   * Fleshes out a path in a table tree. For example, given "fee::foo::fah::koo, this routine will
   * make sure that fee.foo.fah.koo is defined as a table tree structure. The routine will create
   * a table structure where one is needed and leave existing ones alone.
   *
   * At the end of the routine the last table in the structure will be on top of the stack.
   */
  static void FleshOutLuaTableStructure(const std::vector<std::string>& table_names);

  /**
   * Sets a lua constant in the lua state.
   */
  static void SetLuaConstant(const std::string& constant_name, const opensn::Varying& value);

  /**
   * Flushes any commands in the command buffer.
   */
  void FlushConsole();

  /**
   * Generic entry point for wrapper calls.
   */
  static int LuaWrapperCall(lua_State* L);

  /**
   * Makes a formatted output, readible by the documentation scripts, of all the lua wrapper
   * functions.
   */
  void DumpRegister() const;

  /**
   * Given an old status, will update the bindings for only newly registered items.
   */
  void UpdateConsoleBindings(const RegistryStatuses& old_statuses);
};

extern Console& console;

} // namespace opensnlua
