// SPDX-FileCopyrightText: 2024 The OpenSn Authors <https://open-sn.github.io/opensn/>
// SPDX-License-Identifier: MIT

#pragma once

#include "modules/linear_boltzmann_solvers/discrete_ordinates_solver/sweep/spds/spds.h"
#include "modules/linear_boltzmann_solvers/discrete_ordinates_solver/sweep/sweep.h"

namespace opensn
{

class CBC_SPDS : public SPDS
{
public:
  /**
   * Constructs a cell-by-cell sweep-plane data strcture (SPDS) with the given direction and grid.
   *
   * \param omega The angular direction vector.
   * \param grid Reference to the grid.
   * \param allow_cycles Whether cycles are allowed in the local sweep dependency graph.
   */
  CBC_SPDS(const Vector3& omega, const MeshContinuum& grid, bool allow_cycles);

  /// Returns the cell-by-cell task list.
  const std::vector<Task>& GetTaskList() const;

protected:
  /// Cell-by-cell task list.
  std::vector<Task> task_list_;
};

} // namespace opensn
