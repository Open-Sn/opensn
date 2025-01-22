// SPDX-FileCopyrightText: 2024 The OpenSn Authors <https://open-sn.github.io/opensn/>
// SPDX-License-Identifier: MIT

#include "lua/lib/logger.h"
#include "framework/logging/log.h"
#include "framework/runtime.h"

namespace opensnlua
{

void
LuaLogger::Log(opensn::Logger::LOG_LVL level, const std::string& msg)
{
  opensn::log.Log(level) << msg;
}

void
LuaLogger::Log0(const std::string& msg)
{
  opensn::log.Log0() << msg;
}

void
LuaLogger::Log0Warning(const std::string& msg)
{
  opensn::log.Log0Warning() << msg;
}

void
LuaLogger::Log0Error(const std::string& msg)
{
  opensn::log.Log0Error() << msg;
}

void
LuaLogger::Log0Verbose0(const std::string& msg)
{
  opensn::log.Log0Verbose0() << msg;
}

void
LuaLogger::Log0Verbose1(const std::string& msg)
{
  opensn::log.Log0Verbose1() << msg;
}

void
LuaLogger::Log0Verbose2(const std::string& msg)
{
  opensn::log.Log0Verbose2() << msg;
}

void
LuaLogger::LogAll(const std::string& msg)
{
  opensn::log.LogAll() << msg;
}

void
LuaLogger::LogAllWarning(const std::string& msg)
{
  opensn::log.LogAllWarning() << msg;
}

void
LuaLogger::LogAllError(const std::string& msg)
{
  opensn::log.LogAllError() << msg;
}

void
LuaLogger::LogAllVerbose0(const std::string& msg)
{
  opensn::log.LogAllVerbose0() << msg;
}

void
LuaLogger::LogAllVerbose1(const std::string& msg)
{
  opensn::log.LogAllVerbose1() << msg;
}

void
LuaLogger::LogAllVerbose2(const std::string& msg)
{
  opensn::log.LogAllVerbose2() << msg;
}

} // namespace opensnlua
