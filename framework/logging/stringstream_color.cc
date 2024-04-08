#include "framework/logging/stringstream_color.h"
#include "framework/runtime.h"

namespace opensn
{

std::string
StringStreamColor(StringSteamColorCode code)
{
  if (suppress_color)
    return {};
  return std::string("\033[") + std::to_string(code) + "m";
}

} // namespace opensn
