#pragma once

#include "framework/chi_lua.h"

namespace chi_log_utils::lua_utils
{
/** Sets the verbosity level of the Logger.
 * This lua command will overwrite the currently set value.
 *
 * \param int_level int Integer denoting verbosity level. Can be 0,1 or 2 [default:0]
 *
 * \ingroup LuaLogging
 * \author Jan
 */
int chiLogSetVerbosity(lua_State* L);

/**Logs a message depending on the log type specified.
 *
 * \param LogType int Can be any of the log types specified below.
 * \param LogMsg char Message or value to be output to the log.
 *
 * ##_
 *
 * ### LogType
 * LOG_0\n
 * Write a log only if location 0.\n
 *
 * LOG_0WARNING\n
 * Write a log only if location 0 and format it as a warning.\n
 *
 * LOG_0ERROR\n
 * Write a log only if location 0 and format it as an error.\n
 *
 * LOG_0VERBOSE_0\n
 * Same as LOG_0.\n
 *
 * LOG_0VERBOSE_1\n
 * Write a log only if location 0 and the verbosity level is greater
 * or equal to 1.\n
 *
 * LOG_0VERBOSE_2\n
 * Write a log only if location 0 and the verbosity level is greater
 * or equal to 1.\n
 *
 * LOG_ALL, LOG_ALLWARNING, LOG_ALLERROR,\n
 * LOG_ALLVERBOSE_0, LOG_ALLVERBOSE_1, LOG_ALLVERBOSE_2\n
 * Has the same meaning as their LOG_0 counterparts but instead applies to
 * all locations in the parallel context.
 *
 * \ingroup LuaLogging
 * \author Jan
 */
int chiLog(lua_State* L);
/**Processes the sub-events of a repeating event and converts it to a
 * meaningful value (floating-point).
 *
 * \param event_name string Required. Name of the event.
 * \param event_operation_name string Required. What kind of operation to be
 *                            applied. See `chi::ChiLog::EventOperation`
 *
 * \ingroup LuaLogging
 * \return double The processed value.
 */
int chiLogProcessEvent(lua_State* L);
/**Prints the performance graph.
 * \params rank int Optional argument to print the graph for a specific rank.
 *
 * \ingroup LuaLogging
 * */
int chiLogPrintTimingGraph(lua_State* L);
} // namespace chi_log_utils::lua_utils
