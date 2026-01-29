// SPDX-FileCopyrightText: 2025 The OpenSn Authors <https://open-sn.github.io/opensn/>
// SPDX-License-Identifier: MIT

#include "modules/linear_boltzmann_solvers/discrete_ordinates_problem/sweep/angle_set/angle_set.h"
#include <algorithm>

namespace opensn
{

bool
AngleSet::HasAngleIndex(std::uint32_t angle_index) const
{
  return std::find(angles_.begin(), angles_.end(), angle_index) != angles_.end();
}

} // namespace opensn
