// SPDX-FileCopyrightText: 2024 The OpenSn Authors <https://open-sn.github.io/opensn/>
// SPDX-License-Identifier: MIT

#pragma once

#include "modules/linear_boltzmann_solvers/discrete_ordinates_solver/sweep/spds/spds.h"
#include "modules/linear_boltzmann_solvers/discrete_ordinates_solver/sweep/sweep.h"

namespace opensn
{
namespace lbs
{

/**Cell-by-Cell (CBC) Sweep Plane Data Structure*/
class CBC_SPDS : public SPDS
{
public:
  CBC_SPDS(const Vector3& omega,
           const MeshContinuum& grid,
           bool cycle_allowance_flag,
           bool verbose);

  const std::vector<Task>& TaskList() const;

protected:
  std::vector<Task> task_list_;
};

} // namespace lbs
} // namespace opensn
