// SPDX-FileCopyrightText: 2024 The OpenSn Authors <https://open-sn.github.io/opensn/>
// SPDX-License-Identifier: MIT

#include "log.h"
#include "framework/logging/log.h"
#include "framework/console/console.h"
#include "framework/lua.h"
#include "framework/runtime.h"

using namespace opensn;

namespace opensnlua
{
RegisterLuaFunctionInNamespace(LogSetVerbosity, log, SetVerbosity);
RegisterLuaFunctionInNamespace(LogLog, log, Log);

RegisterLuaConstant(LOG_0, Varying(1));
RegisterLuaConstant(LOG_0WARNING, Varying(2));
RegisterLuaConstant(LOG_0ERROR, Varying(3));
RegisterLuaConstant(LOG_0VERBOSE_0, Varying(4));
RegisterLuaConstant(LOG_0VERBOSE_1, Varying(5));
RegisterLuaConstant(LOG_0VERBOSE_2, Varying(6));
RegisterLuaConstant(LOG_ALL, Varying(7));
RegisterLuaConstant(LOG_ALLWARNING, Varying(8));
RegisterLuaConstant(LOG_ALLERROR, Varying(9));
RegisterLuaConstant(LOG_ALLVERBOSE_0, Varying(10));
RegisterLuaConstant(LOG_ALLVERBOSE_1, Varying(11));
RegisterLuaConstant(LOG_ALLVERBOSE_2, Varying(12));

int
LogSetVerbosity(lua_State* L)
{
  const std::string fname = "log.SetVerbosity";
  LuaCheckArgs<int>(L, fname);
  auto level = LuaArg<int>(L, 1);
  if (level <= 2)
  {
    opensn::log.SetVerbosity(level);
  }
  return LuaReturn(L);
}

int
LogLog(lua_State* L)
{
  LuaCheckArgs<int, std::string>(L, "log.Log");
  auto mode = LuaArg<int>(L, 1);
  auto message = LuaArg<std::string>(L, 2);

  opensn::log.Log(static_cast<opensn::Logger::LOG_LVL>(mode)) << message << std::endl;

  return LuaReturn(L);
}

} // namespace opensnlua
