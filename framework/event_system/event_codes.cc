#include "framework/event_system/event_codes.h"
#include <string>
#include <map>

namespace chi
{

Event::EventCode
GetStandardEventCode(const std::string& event_name)
{
  static std::map<std::string, Event::EventCode> event_name_2_code_map{
    {"ProgramStart", Event::ProgramStart},
    {"ProgramExecuted", Event::ProgramExecuted},
    {"SolverPreInitialize", Event::SolverPreInitialize},
    {"SolverInitialized", Event::SolverInitialized},
    {"SolverPreExecution", Event::SolverPreExecution},
    {"SolverExecuted", Event::SolverExecuted},
    {"SolverPreStep", Event::SolverPreStep},
    {"SolverStep", Event::SolverStep},
    {"SolverPreAdvance", Event::SolverPreAdvance},
    {"SolverAdvanced", Event::SolverAdvanced},
  };

  const auto it = event_name_2_code_map.find(event_name);
  if (it != event_name_2_code_map.end()) return it->second;

  return Event::Unknown;
}

} // namespace chi
