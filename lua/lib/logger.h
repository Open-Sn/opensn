// SPDX-FileCopyrightText: 2024 The OpenSn Authors <https://open-sn.github.io/opensn/>
// SPDX-License-Identifier: MIT

#pragma once

#include "framework/logging/log.h"
#include <string>

namespace opensnlua
{

class LuaLogger
{
public:
  static void Log(opensn::Logger::LOG_LVL level, const std::string& msg);
  static void Log0(const std::string& msg);
  static void Log0Warning(const std::string& msg);
  static void Log0Error(const std::string& msg);
  static void Log0Verbose0(const std::string& msg);
  static void Log0Verbose1(const std::string& msg);
  static void Log0Verbose2(const std::string& msg);
  static void LogAll(const std::string& msg);
  static void LogAllWarning(const std::string& msg);
  static void LogAllError(const std::string& msg);
  static void LogAllVerbose0(const std::string& msg);
  static void LogAllVerbose1(const std::string& msg);
  static void LogAllVerbose2(const std::string& msg);
};

extern LuaLogger log;

} // namespace opensnlua
