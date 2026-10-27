// SPDX-FileCopyrightText: 2026 The OpenSn Authors <https://open-sn.github.io/opensn/>
// SPDX-License-Identifier: MIT

#pragma once

#include "modules/linear_boltzmann_solvers/discrete_ordinates_problem/sweep_chunks/gpu_kernel/arguments.h"
#include <cstdint>

namespace opensn::gpu_kernel
{

/**
 * Data for one angle set's ready cells in a fused CBCD sweep kernel launch.
 * A sweep worker combines the ready cells of several of its angle sets into one launch. Each angle
 * set contributes one slot. The slots of a launch occupy consecutive ranges of thread blocks.
 */
struct CBCD_BatchSlot
{
  /// Device pointer to the angle set's sweep arguments.
  const Arguments<SweepKind::CBC>* arguments;
  /// Local IDs of the ready cells, readable by the device.
  const std::uint32_t* cells_to_sweep;
  /// The angle set's saved-angular-flux device buffer (nullptr if fluxes are not saved).
  double* saved_psi;
  /// One past the last thread block of this slot within the launch.
  std::uint64_t block_end;
  /// Number of ready cells.
  std::uint32_t num_cells;
  /// Threads per block along the angle-group axis.
  std::uint32_t block_size_x;
  /// Number of blocks spanning the angle-group axis of one cell row.
  std::uint32_t num_stride_blocks;
};

} // namespace opensn::gpu_kernel
