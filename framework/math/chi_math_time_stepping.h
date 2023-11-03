#pragma once

#include <string>

namespace chi_math
{

enum class SteppingMethod
{
  NONE = 0,
  EXPLICIT_EULER = 1,
  IMPLICIT_EULER = 2,
  CRANK_NICOLSON = 3,
  THETA_SCHEME = 4,
};

std::string SteppingMethodStringName(SteppingMethod);
SteppingMethod SteppingMethodFromString(const std::string& name);

} // namespace chi_math
