// SPDX-FileCopyrightText: 2024 The OpenSn Authors <https://open-sn.github.io/opensn/>
// SPDX-License-Identifier: MIT

#pragma once

#include "framework/parameters/parameter_block.h"
#include "framework/parameters/input_parameters.h"
#include "framework/logging/log_exceptions.h"
#include "framework/utils/utils.h"
extern "C"
{
#include <lua.h>
#include <lualib.h>
#include <lauxlib.h>
}
#include "LuaBridge/LuaBridge.h"
#include <vector>
#include <string>
#include <map>
#include <stack>
#include <functional>

namespace opensnlua
{

/// Class for handling the console and scripting.
class Console
{
  /// Pointer to lua console state
  lua_State* console_state_;
  /// Buffer of commands to execute
  std::vector<std::string> command_buffer_;
  static Console instance_;

  Console() noexcept;

public:
  /// Access to the singleton
  static Console& GetInstance() noexcept;

  lua_State*& GetConsoleState() { return console_state_; }

  std::vector<std::string>& GetCommandBuffer() { return command_buffer_; }

  /// Executes the loop for the console.
  void RunConsoleLoop(char* fileName = nullptr) const;

  /// Executes the given file in the Lua engine.
  int ExecuteFile(const std::string& fileName, int argc, char** argv) const;

  /// Flushes any commands in the command buffer.
  void FlushConsole();

  /**
   * Makes a formatted output, readible by the documentation scripts, of all the lua wrapper
   * functions.
   */
  void DumpRegister() const;

public:
  /**
   * Bind C++ items like classes, functions, etc. to the console
   *
   * \param bind Function with the binding code
   * \return `true`. Need to have a return type so that this can be use with initialization of
   *         static variables
   */
  static bool Bind(std::function<void(lua_State* L)> bind);
};

extern Console& console;

/**
 * Bind a C++ function into a lua namespace.
 *
 * \param namespace_name Lua namespace name
 * \param func_name C++ function name
 */
#define BIND_FUNCTION(namespace_name, func_name)                                                   \
  static bool OpenSnJoinWords(reg_var, __COUNTER__) = opensnlua::Console::Bind(                    \
    [](lua_State* L)                                                                               \
    {                                                                                              \
      luabridge::getGlobalNamespace(L)                                                             \
        .beginNamespace(#namespace_name)                                                           \
        .addFunction(#func_name, func_name)                                                        \
        .endNamespace();                                                                           \
    });

} // namespace opensnlua
