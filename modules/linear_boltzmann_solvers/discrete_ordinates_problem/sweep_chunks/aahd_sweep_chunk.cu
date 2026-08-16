// SPDX-FileCopyrightText: 2025 The OpenSn Authors <https://open-sn.github.io/opensn/>
// SPDX-License-Identifier: MIT

#include "modules/linear_boltzmann_solvers/discrete_ordinates_problem/sweep_chunks/aahd_sweep_chunk.h"
#include "modules/linear_boltzmann_solvers/discrete_ordinates_problem/sweep_chunks/gpu_kernel/main.h"
#include "modules/linear_boltzmann_solvers/discrete_ordinates_problem/sweep_chunks/gpu_kernel/fused_main.h"
#include "modules/linear_boltzmann_solvers/discrete_ordinates_problem/sweep_chunks/gpu_kernel/round_up.h"
#include "modules/linear_boltzmann_solvers/discrete_ordinates_problem/discrete_ordinates_problem.h"
#include "modules/linear_boltzmann_solvers/discrete_ordinates_problem/sweep/angle_set/aahd_angle_set.h"
#include "modules/linear_boltzmann_solvers/discrete_ordinates_problem/sweep/fluds/aahd_fluds.h"
#include "modules/linear_boltzmann_solvers/discrete_ordinates_problem/sweep/spds/aah.h"
#include "modules/linear_boltzmann_solvers/lbs_problem/device/carrier/mesh_carrier.h"
#include "modules/linear_boltzmann_solvers/lbs_problem/device/device_vector_mirror.h"
#include "caliper/cali.h"
#include "caribou/main.hpp"
#include <algorithm>

namespace crb = caribou;

namespace opensn
{

AAHDSweepChunk::AAHDSweepChunk(DiscreteOrdinatesProblem& problem, LBSGroupset& groupset)
  : SweepChunk(problem.GetPhiNewLocal(),
               problem.GetPsiNewLocal()[groupset.id],
               problem.GetGrid(),
               problem.GetSpatialDiscretization(),
               problem.GetUnitCellMatrices(),
               problem.GetCellTransportViews(),
               problem.GetCellOutflowViews(),
               problem.GetQMomentsLocal(),
               groupset,
               problem.GetBlockID2XSMap(),
               problem.GetNumMoments(),
               problem.GetMaxCellDOFCount(),
               problem.GetMinCellDOFCount()),
    problem_(problem)
{
  batch_stream_ = crb::Stream();
}

void
AAHDSweepChunk::Sweep(AngleSet& angle_set)
{
  // prepare arguments
  auto& aahd_angle_set = dynamic_cast<AAHD_AngleSet&>(angle_set);
  auto& fluds = dynamic_cast<AAHD_FLUDS&>(aahd_angle_set.GetFLUDS());
  auto& stream = aahd_angle_set.GetStream();
  gpu_kernel::Arguments<SweepKind::AAH> args(
    problem_, groupset_, aahd_angle_set, fluds, surface_source_active_);
  double* saved_psi = fluds.GetSavedAngularFluxDevicePointer();
  // retrieve SPDS levels
  const auto& spds = dynamic_cast<const AAH_SPDS&>(aahd_angle_set.GetSPDS());
  const auto& levelized_spls = spds.GetLevelizedLocalSubgrid();
  // angle/group count is constant across levels
  auto stride_size = static_cast<unsigned int>(args.flud_data.stride_size);
  for (std::uint32_t level = 0; level < levelized_spls.size(); ++level)
  {
    // Whichever of {level_size, stride_size} is larger drives the fast axis, so warps
    // stay full whichever one is small. Recomputed per level since level_size varies.
    std::size_t level_size = levelized_spls[level].size();
    bool cells_on_fast_axis = level_size >= stride_size;
    unsigned int fast_value = static_cast<unsigned int>(
      std::max<std::size_t>(cells_on_fast_axis ? level_size : stride_size, 1));
    unsigned int slow_value =
      cells_on_fast_axis ? stride_size : static_cast<unsigned int>(level_size);
    unsigned int block_size_x = std::min(gpu_kernel::RoundUp(fast_value), gpu_kernel::threshold);
    unsigned int block_size_y = gpu_kernel::threshold / block_size_x;
    unsigned int grid_size_x = (fast_value + block_size_x - 1) / block_size_x;
    unsigned int grid_size_y = (slow_value + block_size_y - 1) / block_size_y;
    crb::Dim3 block_size(block_size_x, block_size_y);
    crb::Dim3 grid_size(grid_size_x, grid_size_y);
    // perform the sweep on device
    const std::uint32_t* level_data = spds.GetDeviceLevelVector(level);
#if defined(__NVCC__) || defined(__HIPCC__)
    gpu_kernel::SweepKernel<SweepKind::AAH><<<grid_size, block_size, 0, stream>>>(
      args, level_data, static_cast<unsigned int>(level_size), saved_psi, cells_on_fast_axis);
#elif defined(SYCL_LANGUAGE_VERSION) && defined(__INTEL_LLVM_COMPILER)
    stream.parallel_for(
      sycl::nd_range<3>(grid_size * block_size, block_size),
      [=](sycl::nd_item<3> work_index)
      {
        gpu_kernel::SweepKernel<SweepKind::AAH>(
          args, level_data, static_cast<unsigned int>(level_size), saved_psi, cells_on_fast_axis);
      });
#endif
  }
}

void
AAHDSweepChunk::SweepBatchRun(const std::vector<AAHD_AngleSet*>& batch)
{
  if (batch.empty())
    return;

  // Build the shared (angle-set-independent) arguments from any one representative angle set.
  // The kernel replaces its per-angle-set fields with each thread's batch-slot values.
  auto* front_fluds = dynamic_cast<AAHD_FLUDS*>(&batch.front()->GetFLUDS());
  gpu_kernel::Arguments<SweepKind::AAH> shared_args(
    problem_, groupset_, *batch.front(), *front_fluds, surface_source_active_);

  // The run covers the deepest member's level count, and the angle-group axis covers the
  // largest member. Members with fewer levels or angles get empty slots or idle threads; the
  // kernel bounds-checks both against each slot.
  const std::size_t batch_size = batch.size();
  std::size_t max_depth = 0;
  std::size_t max_angleset_size = 0;
  for (auto* as : batch)
  {
    max_depth = std::max(max_depth, as->NumLevels());
    max_angleset_size = std::max(max_angleset_size, as->GetNumAngles());
  }
  if (max_depth == 0)
    return;
  const auto stride_size = static_cast<unsigned int>(max_angleset_size * groupset_.GetNumGroups());

  // Precompute every level's batch-slot data on the host, then copy the whole run to the device
  // at once, so levels launch back-to-back below with no copy or synchronization between them.
  std::vector<gpu_kernel::AAHD_BatchSlot> host_slots(batch_size * max_depth);
  for (std::size_t m = 0; m < batch_size; ++m)
    batch[m]->FillBatchSlots(max_depth, batch_size, host_slots.data() + m);
  std::vector<unsigned int> max_level_size_per_level(max_depth, 0);
  for (std::size_t level = 0; level < max_depth; ++level)
    for (std::size_t m = 0; m < batch_size; ++m)
      max_level_size_per_level[level] =
        std::max(max_level_size_per_level[level], host_slots[level * batch_size + m].num_cells);
  const std::size_t total_slots = host_slots.size();
  if (batch_slots_device_.size() < total_slots)
    batch_slots_device_ = crb::DeviceMemory<gpu_kernel::AAHD_BatchSlot>(total_slots);
  crb::copy(batch_slots_device_, host_slots, total_slots);

  auto* device_slots_base = batch_slots_device_.get();
  for (std::size_t level = 0; level < max_depth; ++level)
  {
    // Same fast-axis reasoning as Sweep() above; grid-z is the batch width.
    unsigned int max_level_size = max_level_size_per_level[level];
    bool cells_on_fast_axis = max_level_size >= stride_size;
    unsigned int fast_value = std::max(cells_on_fast_axis ? max_level_size : stride_size, 1U);
    unsigned int slow_value = std::max(cells_on_fast_axis ? stride_size : max_level_size, 1U);
    unsigned int block_size_x = std::min(gpu_kernel::RoundUp(fast_value), gpu_kernel::threshold);
    unsigned int block_size_y = gpu_kernel::threshold / block_size_x;
    unsigned int grid_size_x = (fast_value + block_size_x - 1) / block_size_x;
    unsigned int grid_size_y = (slow_value + block_size_y - 1) / block_size_y;
    auto grid_size_z = static_cast<unsigned int>(batch_size);
    crb::Dim3 block_size(block_size_x, block_size_y);
    crb::Dim3 grid_size(grid_size_x, grid_size_y, grid_size_z);

    auto* level_slots_ptr = device_slots_base + level * batch_size;
#if defined(__NVCC__) || defined(__HIPCC__)
    gpu_kernel::AAHDFusedSweepKernel<<<grid_size, block_size, 0, batch_stream_>>>(
      shared_args, level_slots_ptr, cells_on_fast_axis);
#elif defined(SYCL_LANGUAGE_VERSION) && defined(__INTEL_LLVM_COMPILER)
    batch_stream_.parallel_for(
      sycl::nd_range<3>(grid_size * block_size, block_size),
      [=](sycl::nd_item<3> work_index)
      { gpu_kernel::AAHDFusedSweepKernel(shared_args, level_slots_ptr, cells_on_fast_axis); });
#endif
  }
  // Exactly one host-blocking wait for the whole run, instead of one per level.
  batch_stream_.synchronize();
}

} // namespace opensn
