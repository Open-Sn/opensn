#include "framework/lua.h"

#include "framework/runtime.h"
#include "framework/logging/log.h"

#include "log_lua.h"
#include "framework/console/console.h"

using namespace opensn;

namespace opensnlua
{
RegisterLuaFunctionAsIs(chiLogSetVerbosity);
RegisterLuaFunctionAsIs(chiLog);
RegisterLuaFunctionAsIs(chiLogProcessEvent);
RegisterLuaFunctionAsIs(chiLogPrintTimingGraph);

RegisterLuaConstantAsIs(LOG_0, Varying(1));
RegisterLuaConstantAsIs(LOG_0WARNING, Varying(2));
RegisterLuaConstantAsIs(LOG_0ERROR, Varying(3));
RegisterLuaConstantAsIs(LOG_0VERBOSE_0, Varying(4));
RegisterLuaConstantAsIs(LOG_0VERBOSE_1, Varying(5));
RegisterLuaConstantAsIs(LOG_0VERBOSE_2, Varying(6));
RegisterLuaConstantAsIs(LOG_ALL, Varying(7));
RegisterLuaConstantAsIs(LOG_ALLWARNING, Varying(8));
RegisterLuaConstantAsIs(LOG_ALLERROR, Varying(9));
RegisterLuaConstantAsIs(LOG_ALLVERBOSE_0, Varying(10));
RegisterLuaConstantAsIs(LOG_ALLVERBOSE_1, Varying(11));
RegisterLuaConstantAsIs(LOG_ALLVERBOSE_2, Varying(12));

int
chiLogSetVerbosity(lua_State* L)
{
  int num_args = lua_gettop(L);

  if (num_args == 0) { return 0; }
  else
  {
    int level = lua_tonumber(L, 1);
    if (level <= 2) { opensn::log.SetVerbosity(level); }
  }
  return 0;
}

int
chiLog(lua_State* L)
{
  int num_args = lua_gettop(L);

  if (num_args != 2) LuaPostArgAmountError("chiLog", 2, num_args);

  int mode = lua_tonumber(L, 1);
  const char* message = lua_tostring(L, 2);

  opensn::log.Log(static_cast<opensn::Logger::LOG_LVL>(mode)) << message << std::endl;

  return 0;
}

int
chiLogProcessEvent(lua_State* L)
{
  const std::string fname = __FUNCTION__;
  const int num_args = lua_gettop(L);
  if (num_args != 2) LuaPostArgAmountError(fname, 2, num_args);

  LuaCheckStringValue(fname, L, 1);
  LuaCheckStringValue(fname, L, 2);

  const std::string event_name = lua_tostring(L, 1);
  const std::string event_operation_name = lua_tostring(L, 2);

  const size_t event_tag = opensn::log.GetExistingRepeatingEventTag(event_name);

  opensn::Logger::EventOperation event_operation;

  if (event_operation_name == "NUMBER_OF_OCCURRENCES")
    event_operation = opensn::Logger::EventOperation::NUMBER_OF_OCCURRENCES;
  else if (event_operation_name == "TOTAL_DURATION")
    event_operation = opensn::Logger::EventOperation::TOTAL_DURATION;
  else if (event_operation_name == "AVERAGE_DURATION")
    event_operation = opensn::Logger::EventOperation::AVERAGE_DURATION;
  else if (event_operation_name == "MAX_VALUE")
    event_operation = opensn::Logger::EventOperation::MAX_VALUE;
  else if (event_operation_name == "AVERAGE_VALUE")
    event_operation = opensn::Logger::EventOperation::AVERAGE_VALUE;
  else
    ChiInvalidArgument("Unsupported event operation name \"" + event_operation_name + "\".");

  const double value = opensn::log.ProcessEvent(event_tag, event_operation);

  lua_pushnumber(L, static_cast<lua_Number>(value));
  return 1;
}

int
chiLogPrintTimingGraph(lua_State* L)
{
  const std::string fname = __FUNCTION__;
  const int num_args = lua_gettop(L);
  auto& chitech_timing = opensn::log.GetTimingBlock("ChiTech");

  int rank = 0;
  if (num_args >= 1)
  {
    LuaCheckIntegerValue(fname, L, 1);
    rank = lua_tointeger(L, 1);
  }

  ChiInvalidArgumentIf(rank >= opensn::mpi.process_count,
                       "rank >= process_count, i.e., " + std::to_string(rank) +
                         " >= " + std::to_string(opensn::mpi.process_count));

  if (opensn::mpi.location_id == rank)
    opensn::log.LogAll() << "\nPerformance Graph:\n" << chitech_timing.MakeGraphString();

  return 0;
}

} // namespace opensnlua
