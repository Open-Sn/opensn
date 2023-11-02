#include <string>
#include <map>

namespace chi
{

int GetStandardEventCode(const std::string& event_name)
{
  static std::map<std::string, int> event_name_2_code_map{
    {"ProgramStart", 1},
    {"ProgramExecuted", 2},
    {"SolverPreInitialize", 31},
    {"SolverInitialized", 32},
    {"SolverPreExecution", 33},
    {"SolverExecuted", 34},
    {"SolverPreStep", 35},
    {"SolverStep", 36},
    {"SolverPreAdvance", 37},
    {"SolverAdvanced", 38},
  };

  const auto it = event_name_2_code_map.find(event_name);
  if (it != event_name_2_code_map.end()) return it->second;

  return 0;
}

} // namespace chi