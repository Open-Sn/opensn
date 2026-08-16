// SPDX-FileCopyrightText: 2026 The OpenSn Authors <https://open-sn.github.io/opensn/>
// SPDX-License-Identifier: MIT

#pragma once

#include "modules/linear_boltzmann_solvers/discrete_ordinates_problem/sweep/fluds/aahd_structs.h"
#include <cstdint>

namespace opensn::gpu_kernel
{

// Per-angle-set data for one grid-z slot of a fused (multi-angle-set) AAHD sweep kernel launch.
// Fields shared by every angle set in a groupset's sweep (mesh/quadrature/source/phi/boundary/
// groupset info) are passed once per launch via `Arguments<SweepKind::AAH>`. Only the fields
// that vary per angle set live here.
struct AAHD_BatchSlot
{
  /// Device pointer to this angle set's direction index (indices, for angleset_size > 1).
  const std::uint32_t* __restrict__ directions = nullptr;
  /// Device pointer to this angle set's boundary offset table.
  const std::uint64_t* __restrict__ boundary_offset = nullptr;
  /// Number of angles in this angle set.
  std::uint32_t angleset_size = 0;
  /// Device pointer to this angle set's FLUDS common node-index table.
  const std::uint64_t* __restrict__ flud_index = nullptr;
  /// This angle set's device FLUDS pointer set.
  AAHD_FLUDSPointerSet flud_data;
  /// Device pointer to this level's cell list (nullptr for an empty slot).
  const std::uint32_t* __restrict__ cells_to_sweep = nullptr;
  /// Number of cells in this level (0 for an empty slot).
  unsigned int num_cells = 0;
  /// This angle set's saved-angular-flux device buffer (nullptr if not in use).
  double* __restrict__ saved_psi = nullptr;
};

} // namespace opensn::gpu_kernel
