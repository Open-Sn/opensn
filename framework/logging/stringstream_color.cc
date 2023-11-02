#include "stringstream_color.h"

#include "chi_runtime.h"

std::string
chi::StringStreamColor(StringSteamColorCode code)
{
  if (Chi::run_time::suppress_color_) return {};
  return std::string("\033[") + std::to_string(code) + "m";
}