// SPDX-FileCopyrightText: 2024 The OpenSn Authors <https://open-sn.github.io/opensn/>
// SPDX-License-Identifier: MIT

#pragma once

#include "modules/linear_boltzmann_solvers/discrete_ordinates_problem/sweep/spds/spds.h"
#include "modules/linear_boltzmann_solvers/discrete_ordinates_problem/sweep/sweep.h"

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
  CBC_SPDS(const Vector3& omega, const std::shared_ptr<MeshContinuum>& grid, bool allow_cycles);

  /// Returns the cell-by-cell task list.
  const std::vector<Task>& GetTaskList() const;

  /// Returns the minimum number of slots needed for the pool allocator in CBC_FLUDS.
  size_t GetMinNumPoolAllocatorSlots() const { return min_num_pool_allocator_slots_; }

private:
  /**
   * Simulates a sweep over the local cells to calculate the minimum number of slots needed for
   * the pool allocator in CBC_FLUDS.
   * Due to the asynchronous nature of communication of the CBC algorithm, the simulated sweep sets
   * aside a slot for each cell that has either remote upwind or downwind dependencies, which cannot
   * be reused during a sweep. For cells that have only local dependencies, slots can be reused.
   */
  size_t SimulateLocalSweep() const;

protected:
  /// Cell-by-cell task list.
  std::vector<Task> task_list_;

  /// Minimum number of slots needed for CBC_FLUDS pool allocator.
  size_t min_num_pool_allocator_slots_ = 0;
};

} // namespace opensn
