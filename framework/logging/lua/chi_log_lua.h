#ifndef CHITECH_CHI_LOG_LUA_H
#define CHITECH_CHI_LOG_LUA_H

#include "chi_lua.h"

namespace chi_log_utils::lua_utils
{
int chiLogSetVerbosity(lua_State* L);
int chiLog(lua_State* L);
int chiLogProcessEvent(lua_State* L);
int chiLogPrintTimingGraph(lua_State* L);
} // namespace chi_log_utils::lua_utils

#endif // CHITECH_CHI_LOG_LUA_H
