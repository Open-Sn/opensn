// SPDX-FileCopyrightText: 2024 The OpenSn Authors <https://open-sn.github.io/opensn/>
// SPDX-License-Identifier: MIT

#pragma once

#include "modules/linear_boltzmann_solvers/discrete_ordinates_problem/sweep/angle_set/angle_set.h"

namespace opensn
{

/// Manages the workstages of a single angleset group.
class AngleSetGroup
{
public:
  std::vector<std::shared_ptr<AngleSet>>& GetAngleSets() { return angle_sets_; }

private:
  std::vector<std::shared_ptr<AngleSet>> angle_sets_;
};

} // namespace opensn
