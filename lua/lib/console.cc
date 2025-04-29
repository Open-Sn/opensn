// SPDX-FileCopyrightText: 2024 The OpenSn Authors <https://open-sn.github.io/opensn/>
// SPDX-License-Identifier: MIT

#include "lua/lib/console.h"
#include "framework/runtime.h"
#include "framework/logging/log.h"
#include "framework/logging/log_exceptions.h"
#include "framework/utils/utils.h"
#include "framework/runtime.h"
#include "config.h"
#include <iostream>

using namespace opensn;

namespace opensnlua
{

Console& console = Console::GetInstance();

Console&
Console::GetInstance() noexcept
{
  static Console singleton;
  return singleton;
}

Console::Console() noexcept : console_state_(luaL_newstate())
{
  luaL_openlibs(console_state_);
}

void
Console::FlushConsole()
{
  try
  {
    for (auto& command : command_buffer_)
    {
      bool error = luaL_dostring(console_state_, command.c_str());
      if (error)
      {
        opensn::log.LogAll() << lua_tostring(console_state_, -1);
        lua_pop(console_state_, 1);
      }
    }
  }
  catch (const std::exception& e)
  {
    opensn::log.LogAllError() << e.what();
    opensn::mpi_comm.abort(EXIT_FAILURE);
  }
}

void
Console::RunConsoleLoop(char*) const
{
  opensn::log.Log() << "Console loop started. "
                    << "Type \"exit\" to quit (or Ctl-C).";

  /** Executes a string within the lua-console. */
  auto LuaDoString = [this](const std::string& the_string)
  {
    bool error = luaL_dostring(console_state_, the_string.c_str());
    if (error)
    {
      opensn::log.LogAll() << lua_tostring(console_state_, -1);
      lua_pop(console_state_, 1);
    }
  };

  const bool HOME = opensn::mpi_comm.rank() == 0;

  while (true)
  {
    std::string console_input;

    if (HOME)
      std::cin >> console_input; // Home will be waiting here

    mpi_comm.broadcast(console_input, 0);
    if (console_input == "exit")
      break;

    try
    {
      LuaDoString(console_input);
    }
    catch (const std::exception& e)
    {
      opensn::log.LogAllError() << e.what();
      opensn::mpi_comm.abort(EXIT_FAILURE);
    }
  } // while not termination posted

  opensn::log.Log() << "Console loop stopped successfully.";
}

int
Console::ExecuteFile(const std::string& fileName, int argc, char** argv) const
{
  lua_State* L = this->console_state_;
  if (not fileName.empty())
  {
    if (argc > 0)
    {
      std::vector<std::string> args;
      args.resize(argc);
      for (int i = 0; i < argc; ++i)
        args[i] = std::string(argv[i]);
    }
    int error = luaL_dofile(this->console_state_, fileName.c_str());

    if (error > 0)
    {
      opensn::log.LogAllError() << lua_tostring(this->console_state_, -1);
      return EXIT_FAILURE;
    }
  }
  return EXIT_SUCCESS;
}

void
Console::DumpRegister() const
{
}

bool
Console::Bind(std::function<void(lua_State* L)> bind)
{
  bind(GetInstance().GetConsoleState());
  return true;
}

} // namespace opensnlua
