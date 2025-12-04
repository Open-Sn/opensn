// SPDX-FileCopyrightText: 2025 The OpenSn Authors <https://open-sn.github.io/opensn/>
// SPDX-License-Identifier: MIT

#include "modules/linear_boltzmann_solvers/discrete_ordinates_problem/sweep/angle_set/cbc_angle_set.h"
#include "modules/linear_boltzmann_solvers/discrete_ordinates_problem/sweep/angle_set/cbc_angle_set_helpers.h"
#include "caribou/cuda/stream.hpp"

namespace opensn
{

caribou::Stream&
GetCBCAngleSetStream(CBC_AngleSet& angle_set)
{
  return *static_cast<caribou::Stream*>(angle_set.GetStream());
}

const caribou::Stream&
GetCBCAngleSetStream(const CBC_AngleSet& angle_set)
{
  return *static_cast<const caribou::Stream*>(angle_set.GetStream());
}

} // namespace opensn