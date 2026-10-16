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

  // Fused-kernel launch parameters in units of the device's warp (wavefront) size and
  // compute-unit count, so they scale across NVIDIA, AMD, and Intel GPUs.
  const std::size_t warp_size = crb::get_warp_size();
#if defined(__NVCC__) || defined(__HIPCC__)
  const auto num_compute_units = static_cast<std::size_t>(crb::impl::get_device_attribute(
    crb::impl::get_current_device(),
    ::CUDA_OR_HIP(DevAttrMultiProcessorCount, DeviceAttributeMultiprocessorCount)));
#elif defined(SYCL_LANGUAGE_VERSION) && defined(__INTEL_LLVM_COMPILER)
  const auto num_compute_units = static_cast<std::size_t>(
    batch_stream_.get_device().get_info<::sycl::info::device::max_compute_units>());
#else
  const std::size_t num_compute_units = 1;
#endif
  fused_launch_.linear_block_size = static_cast<unsigned int>(2 * warp_size);
  fused_launch_.min_small_level_stride = static_cast<unsigned int>(warp_size / 4);
  fused_launch_.linear_mapping_min_threads = 4 * warp_size * num_compute_units;
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
  // compute block size
  unsigned int stride_size =
    gpu_kernel::RoundUp(static_cast<unsigned int>(args.flud_data.stride_size));
  unsigned int block_size_x = std::min(stride_size, gpu_kernel::threshold);
  unsigned int block_size_y = gpu_kernel::threshold / block_size_x;
  crb::Dim3 block_size(block_size_x, block_size_y);
  unsigned int grid_size_x = (stride_size + gpu_kernel::threshold - 1) / gpu_kernel::threshold;
  for (std::uint32_t level = 0; level < levelized_spls.size(); ++level)
  {
    // compute grid size
    std::size_t level_size = levelized_spls[level].size();
    unsigned int grid_size_y = (level_size + block_size_y - 1) / block_size_y;
    crb::Dim3 grid_size(grid_size_x, grid_size_y);
    // perform the sweep on device
    const std::uint32_t* level_data = spds.GetDeviceLevelVector(level);
#if defined(__NVCC__) || defined(__HIPCC__)
    gpu_kernel::SweepKernel<SweepKind::AAH><<<grid_size, block_size, 0, stream>>>(
      args, level_data, static_cast<unsigned int>(level_size), saved_psi);
#elif defined(SYCL_LANGUAGE_VERSION) && defined(__INTEL_LLVM_COMPILER)
    stream.parallel_for(sycl::nd_range<3>(grid_size * block_size, block_size),
                        [=](sycl::nd_item<3> work_index)
                        {
                          gpu_kernel::SweepKernel<SweepKind::AAH>(
                            args, level_data, static_cast<unsigned int>(level_size), saved_psi);
                        });
#endif
  }
}

namespace
{

/// Launches one level of a fused run on `stream`.
template <std::uint32_t fixed_ndofs>
void
LaunchFusedLevel(crb::Stream& stream,
                 const gpu_kernel::Arguments<SweepKind::AAH>& shared_args,
                 const gpu_kernel::AAHD_BatchSlot* level_slots,
                 crb::Dim3 grid,
                 crb::Dim3 block,
                 bool linear_mapping)
{
#if defined(__NVCC__) || defined(__HIPCC__)
  gpu_kernel::AAHDFusedSweepKernel<fixed_ndofs>
    <<<grid, block, 0, stream>>>(shared_args, level_slots, linear_mapping);
#elif defined(SYCL_LANGUAGE_VERSION) && defined(__INTEL_LLVM_COMPILER)
  stream.parallel_for(
    sycl::nd_range<3>(grid * block, block),
    [=](sycl::nd_item<3> work_index)
    { gpu_kernel::AAHDFusedSweepKernel<fixed_ndofs>(shared_args, level_slots, linear_mapping); });
#endif
}

} // namespace

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

  // Compile for the cell DOF count when every local cell has the same one.
  const std::uint32_t fixed_ndofs =
    (min_num_cell_dofs_ == max_num_cell_dofs_) ? max_num_cell_dofs_ : 0;

  auto* device_slots_base = batch_slots_device_.get();
  for (std::size_t level = 0; level < max_depth; ++level)
  {
    // Large levels use the linear (cell, angle-group) mapping, which leaves no idle lanes. Small
    // levels put the angle-groups on x and the cells on y, which spreads their few threads over
    // more warps for latency hiding, unless a cell has too few angle-groups to fill a useful part
    // of a warp. grid-z is the batch width.
    const unsigned int max_level_size = std::max(max_level_size_per_level[level], 1U);
    const std::size_t level_threads =
      static_cast<std::size_t>(max_level_size) * stride_size * batch_size;
    const bool linear_mapping = level_threads >= fused_launch_.linear_mapping_min_threads or
                                stride_size < fused_launch_.min_small_level_stride;
    crb::Dim3 block_size;
    crb::Dim3 grid_size;
    if (linear_mapping)
    {
      const unsigned int linear_block_size = fused_launch_.linear_block_size;
      block_size = crb::Dim3(linear_block_size, 1);
      grid_size =
        crb::Dim3((max_level_size * stride_size + linear_block_size - 1) / linear_block_size,
                  1,
                  static_cast<unsigned int>(batch_size));
    }
    else
    {
      unsigned int block_size_x = std::min(gpu_kernel::RoundUp(stride_size), gpu_kernel::threshold);
      unsigned int block_size_y = gpu_kernel::threshold / block_size_x;
      block_size = crb::Dim3(block_size_x, block_size_y);
      grid_size = crb::Dim3((stride_size + block_size_x - 1) / block_size_x,
                            (max_level_size + block_size_y - 1) / block_size_y,
                            static_cast<unsigned int>(batch_size));
    }

    auto* level_slots_ptr = device_slots_base + level * batch_size;
    auto launch = [&](auto fixed_c)
    {
      LaunchFusedLevel<decltype(fixed_c)::value>(
        batch_stream_, shared_args, level_slots_ptr, grid_size, block_size, linear_mapping);
    };
    switch (fixed_ndofs)
    {
      case 2:
        launch(std::integral_constant<std::uint32_t, 2>{});
        break;
      case 3:
        launch(std::integral_constant<std::uint32_t, 3>{});
        break;
      case 4:
        launch(std::integral_constant<std::uint32_t, 4>{});
        break;
      case 6:
        launch(std::integral_constant<std::uint32_t, 6>{});
        break;
      case 8:
        launch(std::integral_constant<std::uint32_t, 8>{});
        break;
      default:
        launch(std::integral_constant<std::uint32_t, 0>{});
        break;
    }
  }
  // Exactly one host-blocking wait for the whole run, instead of one per level.
  batch_stream_.synchronize();
}

} // namespace opensn
