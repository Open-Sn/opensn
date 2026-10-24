// SPDX-FileCopyrightText: 2026 The OpenSn Authors <https://open-sn.github.io/opensn/>
// SPDX-License-Identifier: MIT

#include "modules/linear_boltzmann_solvers/discrete_ordinates_problem/sweep_chunks/cbcd_sweep_chunk.h"
#include "modules/linear_boltzmann_solvers/discrete_ordinates_problem/sweep_chunks/gpu_kernel/cbcd_fused_main.h"
#include "modules/linear_boltzmann_solvers/discrete_ordinates_problem/sweep_chunks/gpu_kernel/round_up.h"
#include "caliper/cali.h"
#include <algorithm>
#include <cstdint>
#include <vector>

namespace opensn
{

CBCDSweepChunk::CBCDSweepChunk(DiscreteOrdinatesProblem& problem, LBSGroupset& groupset)
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
  const auto num_angle_sets = groupset.angle_agg->GetNumAngleSets();
  angle_sets_.reserve(num_angle_sets);
  angle_set_launches_.reserve(num_angle_sets);
  for (auto& as : *(groupset.angle_agg))
  {
    auto* angle_set = dynamic_cast<CBCD_AngleSet*>(as.get());
    auto* fluds = dynamic_cast<CBCD_FLUDS*>(&(angle_set->GetFLUDS()));
    angle_sets_.push_back(angle_set);

    fluds->SetSaveAngularFlux(problem_.SaveAngularFluxEnabled());
    angle_set->RefreshDeviceData();

    // This is the same block shape as AAHDSweepChunk::Sweep, with x covering the angle-group stride
    // and y covering the cells.
    const auto stride_size = gpu_kernel::RoundUp(static_cast<unsigned int>(fluds->GetStrideSize()));
    const auto block_size_x = std::min(stride_size, gpu_kernel::threshold);
    const auto block_size_y = gpu_kernel::threshold / block_size_x;
    const auto grid_size_x = (stride_size + gpu_kernel::threshold - 1) / gpu_kernel::threshold;
    angle_set_launches_.push_back({crb::Dim3(block_size_x, block_size_y),
                                   grid_size_x,
                                   fluds,
                                   fluds->GetSavedAngularFluxDevicePointer()});
  }
  // RefreshKernelArguments builds the arguments before each sweep.
  host_arguments_.reserve(angle_sets_.size());

  const auto& boundaries = problem_.GetSweepBoundaries();
  needs_boundary_upload_ =
    std::any_of(boundaries.begin(),
                boundaries.end(),
                [&groupset](const auto& entry)
                { return entry.second->CountDelayedAngularDOFsOld(groupset.id) != 0; });

  if (not angle_sets_.empty())
  {
    std::vector<AngleSet*> base_angle_sets(angle_sets_.begin(), angle_sets_.end());
    async_comm_ = std::make_unique<CBCD_AsynchronousCommunicator>(
      base_angle_sets, angle_sets_.front()->GetSweepCommunicator());
    for (auto* angle_set : angle_sets_)
      angle_set->SetCommunicator(*async_comm_);
  }
}

CBCDSweepChunk::~CBCDSweepChunk()
{
  StopCommunicator();
}

void
CBCDSweepChunk::StartCommunicator(const std::size_t num_workers)
{
  if (lanes_.size() != num_workers)
  {
    // A launch has at most one batch slot per angle set of its worker.
    std::vector<std::size_t> num_owned(num_workers, 0);
    for (std::size_t angle_set_id = 0; angle_set_id < angle_sets_.size(); ++angle_set_id)
      ++num_owned[CBCD_AngleSetWorker(angle_set_id, num_workers)];
    const auto max_owned = *std::max_element(num_owned.begin(), num_owned.end());
    lanes_.clear();
    for (std::size_t worker = 0; worker < num_workers; ++worker)
    {
      auto lane = std::make_unique<FusedLane>();
      lane->batch_slots.reserve(max_owned);
      lanes_.push_back(std::move(lane));
    }
  }
  if (async_comm_)
    async_comm_->Start(num_workers);
}

void
CBCDSweepChunk::StopCommunicator()
{
  if (async_comm_)
    async_comm_->Stop();
}

void
CBCDSweepChunk::RefreshKernelArguments()
{
  CALI_CXX_MARK_SCOPE("CBCDSweepChunk::RefreshKernelArguments");

  host_arguments_.clear();
  for (std::size_t angle_set_id = 0; angle_set_id < angle_sets_.size(); ++angle_set_id)
  {
    auto& launch = angle_set_launches_[angle_set_id];
    host_arguments_.push_back(gpu_kernel::Arguments<SweepKind::CBC>(
      problem_, groupset_, *angle_sets_[angle_set_id], *launch.fluds, IsSurfaceSourceActive()));
    host_arguments_.back().flud_data.nonlocal_incoming_psi = async_comm_->GetIncomingPacketValues();
    launch.device_saved_psi = launch.fluds->GetSavedAngularFluxDevicePointer();
  }
  if (device_arguments_.size() != host_arguments_.size())
    device_arguments_ =
      crb::DeviceMemory<gpu_kernel::Arguments<SweepKind::CBC>>(host_arguments_.size());
  crb::copy(device_arguments_, host_arguments_, host_arguments_.size());
}

void
CBCDSweepChunk::AddBatchSlot(const std::size_t worker_id,
                             const std::size_t angle_set_id,
                             const std::uint32_t* local_cell_ids,
                             const std::uint32_t num_ready_cells)
{
  CALI_CXX_MARK_SCOPE("CBCDSweepChunk::AddBatchSlot");

  auto& lane = *lanes_[worker_id];
  auto& launch = angle_set_launches_[angle_set_id];
  const auto num_cell_blocks =
    (std::uint64_t{num_ready_cells} + launch.threads_per_block.y - 1) / launch.threads_per_block.y;
  const auto block_begin = lane.batch_slots.empty() ? 0 : lane.batch_slots.back().block_end;
  lane.batch_slots.push_back({device_arguments_.get() + angle_set_id,
                              local_cell_ids,
                              launch.device_saved_psi,
                              block_begin + num_cell_blocks * launch.num_stride_blocks,
                              num_ready_cells,
                              launch.threads_per_block.x,
                              launch.num_stride_blocks});
}

void
CBCDSweepChunk::SubmitFusedLaunch(const std::size_t worker_id)
{
  auto& lane = *lanes_[worker_id];
  if (lane.in_flight or lane.batch_slots.empty())
    return;
  CALI_CXX_MARK_SCOPE("CBCDSweepChunk::SubmitFusedLaunch");
  const auto num_blocks = lane.batch_slots.back().block_end;
  const auto num_batch_slots = static_cast<std::uint32_t>(lane.batch_slots.size());
  const auto* batch_slots = lane.batch_slots.data();
  // The grid can be capped at 65535 blocks because the kernel strides over the blocks.
  crb::Dim3 grid_size(static_cast<unsigned int>(std::min<std::uint64_t>(num_blocks, 65535)));
  crb::Dim3 block_size(gpu_kernel::threshold);
  auto& stream = lane.stream;
#if defined(__NVCC__) || defined(__HIPCC__)
  gpu_kernel::CBCDFusedSweepKernel<<<grid_size, block_size, 0, stream>>>(
    batch_slots, num_batch_slots, num_blocks);
#elif defined(SYCL_LANGUAGE_VERSION) && defined(__INTEL_LLVM_COMPILER)
  stream.parallel_for(
    sycl::nd_range<3>(grid_size * block_size, block_size),
    [=](sycl::nd_item<3> work_index)
    { gpu_kernel::CBCDFusedSweepKernel(batch_slots, num_batch_slots, num_blocks); });
#endif
  lane.in_flight = true;
}

} // namespace opensn
