// SPDX-FileCopyrightText: 2026 The OpenSn Authors <https://open-sn.github.io/opensn/>
// SPDX-License-Identifier: MIT

#pragma once

#include "modules/linear_boltzmann_solvers/discrete_ordinates_problem/sweep_chunks/gpu_kernel/aahd_batch_slot.h"
#include "modules/linear_boltzmann_solvers/discrete_ordinates_problem/sweep_chunks/gpu_kernel/arguments.h"
#include "modules/linear_boltzmann_solvers/discrete_ordinates_problem/sweep_chunks/gpu_kernel/main.h"
#include "caribou/main.hpp"

namespace crb = caribou;

namespace opensn::gpu_kernel
{

// Fused AAHD sweep kernel. Each grid-z slice (blockIdx.z / SYCL dimension 0) is a different
// angle set's current level, so one launch covers many angle sets' worth of work instead of
// one. Separate from `SweepKernel` in main.h (used by the default, non-fused scheduler); the
// per-thread compute below is the same `SweepDispatch`/solver.h logic, just fed an
// `Arguments<AAH>` assembled from the shared launch data plus this thread's batch slot.
__CRB_GLOBAL_FUNC__ void
AAHDFusedSweepKernel(Arguments<SweepKind::AAH> shared_args,
                     const AAHD_BatchSlot* batch_slots,
                     bool cells_on_fast_axis)
{
#if defined(__NVCC__) || defined(__HIPCC__)
  unsigned int slot_idx = blockIdx.z;
  unsigned int fast_idx = threadIdx.x + blockDim.x * blockIdx.x;
  unsigned int slow_idx = threadIdx.y + blockDim.y * blockIdx.y;
#elif defined(SYCL_LANGUAGE_VERSION) && defined(__INTEL_LLVM_COMPILER)
  auto work_index = ::sycl::ext::oneapi::this_work_item::get_nd_item<3>();
  unsigned int slot_idx = work_index.get_global_id(0);
  unsigned int fast_idx = work_index.get_global_id(2);
  unsigned int slow_idx = work_index.get_global_id(1);
#endif
  const AAHD_BatchSlot& slot = batch_slots[slot_idx];
  unsigned int cell_idx = cells_on_fast_axis ? fast_idx : slow_idx;
  unsigned int angle_group_idx = cells_on_fast_axis ? slow_idx : fast_idx;
  if (cell_idx >= slot.num_cells ||
      angle_group_idx >= shared_args.groupset_size * slot.angleset_size)
    return;

  // Assemble this thread's Arguments<AAH> from the shared launch data plus this slot's
  // per-angle-set fields.
  Arguments<SweepKind::AAH> args = shared_args;
  args.directions = slot.directions;
  args.boundary_offset = slot.boundary_offset;
  args.angleset_size = slot.angleset_size;
  args.flud_index = slot.flud_index;
  args.flud_data = slot.flud_data;

  unsigned int angle_idx = angle_group_idx / args.groupset_size;
  unsigned int group_idx = angle_group_idx - angle_idx * args.groupset_size;
  const std::uint32_t cell_local_idx = slot.cells_to_sweep[cell_idx];
  CellView cell;
  MeshView(args.mesh_data).GetCellView(cell, cell_local_idx);
  if (cell.num_nodes == 0)
    return;
  auto [cell_edge_data, _] = GetCellDataIndex(args.flud_index, cell_local_idx);
  std::uint32_t num_moments;
  std::uint32_t direction_num = args.directions[angle_idx];
  DirectionView direction;
  {
    QuadratureView quadrature(args.quad_data);
    num_moments = quadrature.num_moments;
    quadrature.GetDirectionView(direction, direction_num);
  }
  opensn::gpu_kernel::SweepDispatch<SweepKind::AAH>(cell.num_nodes,
                                                    args,
                                                    cell,
                                                    direction,
                                                    cell_edge_data,
                                                    angle_group_idx,
                                                    group_idx,
                                                    num_moments,
                                                    slot.saved_psi);
}

} // namespace opensn::gpu_kernel
