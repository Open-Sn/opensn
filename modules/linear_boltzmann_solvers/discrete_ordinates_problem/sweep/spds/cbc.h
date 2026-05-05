// SPDX-FileCopyrightText: 2024 The OpenSn Authors <https://open-sn.github.io/opensn/>
// SPDX-License-Identifier: MIT

#pragma once

#include "modules/linear_boltzmann_solvers/discrete_ordinates_problem/sweep/spds/spds.h"
#include "modules/linear_boltzmann_solvers/discrete_ordinates_problem/sweep/sweep.h"

namespace opensn
{

/// CBC sweep-plane data structure.
class CBC_SPDS : public SPDS
{
public:
  CBC_SPDS(const Vector3& omega, const std::shared_ptr<MeshContinuum>& grid, bool allow_cycles);

  const std::vector<Task>& GetTaskList() const;

protected:
  /// Cell-by-cell task list.
  std::vector<Task> task_list_;
};

} // namespace opensn
