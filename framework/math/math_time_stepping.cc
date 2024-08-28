// SPDX-FileCopyrightText: 2024 The OpenSn Authors <https://open-sn.github.io/opensn/>
// SPDX-License-Identifier: MIT

#include "framework/math/math_time_stepping.h"
#include <stdexcept>

namespace opensn
{

#pragma GCC diagnostic push
#pragma GCC diagnostic error "-Wswitch-enum"

std::string
SteppingMethodStringName(SteppingMethod method)
{
  switch (method)
  {
    case SteppingMethod::NONE:
      return "none";
    case SteppingMethod::EXPLICIT_EULER:
      return "explicit_euler";
    case SteppingMethod::IMPLICIT_EULER:
      return "implicit_euler";
    case SteppingMethod::CRANK_NICOLSON:
      return "crank_nicholson";
    case SteppingMethod::THETA_SCHEME:
      return "theta_scheme";
    default:
      throw std::logic_error(__PRETTY_FUNCTION__);
  }
}
#pragma GCC diagnostic pop

} // namespace opensn
