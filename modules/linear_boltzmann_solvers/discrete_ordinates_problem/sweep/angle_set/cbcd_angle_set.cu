// SPDX-FileCopyrightText: 2026 The OpenSn Authors <https://open-sn.github.io/opensn/>
// SPDX-License-Identifier: MIT

#include "modules/linear_boltzmann_solvers/discrete_ordinates_problem/sweep/angle_set/cbcd_angle_set.h"
#include "modules/linear_boltzmann_solvers/discrete_ordinates_problem/sweep/communicators/cbcd_async_comm.h"
#include "modules/linear_boltzmann_solvers/discrete_ordinates_problem/sweep/fluds/cbcd_fluds.h"
#include "modules/linear_boltzmann_solvers/discrete_ordinates_problem/sweep/spds/cbc.h"
#include "modules/linear_boltzmann_solvers/discrete_ordinates_problem/sweep_chunks/cbcd_sweep_chunk.h"
#include "modules/linear_boltzmann_solvers/discrete_ordinates_problem/sweep/boundary/sweep_boundary.h"
#include "modules/linear_boltzmann_solvers/lbs_problem/groupset/lbs_groupset.h"
#include "framework/mesh/mesh/mesh.h"
#include "framework/utils/error.h"
#include "caliper/cali.h"
#include <algorithm>

namespace opensn
{

CBCD_AngleSet::CBCD_AngleSet(size_t id,
                             const LBSGroupset& groupset,
                             const SPDS& spds,
                             std::shared_ptr<FLUDS>& fluds,
                             const std::vector<size_t>& angle_indices,
                             std::map<uint64_t, std::shared_ptr<SweepBoundary>>& boundaries,
                             const SweepCommunicator& sweep_communicator)
  : AngleSet(id, groupset, spds, fluds, angle_indices, boundaries),
    cbc_spds_(dynamic_cast<const CBC_SPDS&>(spds)),
    sweep_communicator_(sweep_communicator),
    cbcd_fluds_(dynamic_cast<CBCD_FLUDS&>(*fluds_)),
    stream_(),
    device_angle_indices_(angles_.size())
{
  cbcd_fluds_.GetStream() = stream_;
  for (auto& buffer : batch_pipeline_.cells)
    buffer.reserve(cbc_spds_.GetTaskList().size());
  BuildCellTaskGraph();
}

void
CBCD_AngleSet::SyncDeviceAngleIndices()
{
  crb::copy(device_angle_indices_, angles_, angles_.size());
}

void
CBCD_AngleSet::RefreshDeviceData()
{
  SyncDeviceAngleIndices();
  BuildReflectingCellMask();
  cbcd_fluds_.BuildOutgoingBoundaryMap(boundaries_, *this);
  cbcd_fluds_.BuildIncomingBoundaryMap(boundaries_, *this);
}

void
CBCD_AngleSet::ResetSweepDependencies()
{
  unresolved_sweep_dependencies_.store(num_dependencies_, std::memory_order_relaxed);
}

void
CBCD_AngleSet::BuildReflectingCellMask()
{
  const auto& common_data = cbcd_fluds_.GetCommonData();
  const auto num_cells = common_data.GetNumLocalCells();
  has_outgoing_reflecting_face_.assign(num_cells, 0);
  num_reflecting_cells_ = 0;
  for (std::size_t cell_local_id = 0; cell_local_id < num_cells; ++cell_local_id)
  {
    for (const auto& face : common_data.GetOutgoingBoundaryFaces(cell_local_id))
    {
      const auto boundary_it = boundaries_.find(face.boundary_id);
      if (boundary_it != boundaries_.end() and boundary_it->second->IsReflecting())
      {
        has_outgoing_reflecting_face_[cell_local_id] = 1;
        ++num_reflecting_cells_;
        break;
      }
    }
  }
}

void
CBCD_AngleSet::BuildCellTaskGraph()
{
  const auto& task_list = cbc_spds_.GetTaskList();
  const auto num_cells = task_list.size();

  remaining_cell_dependencies_.resize(num_cells);
  cell_successor_offsets_.assign(num_cells + 1, 0);
  for (std::size_t task_idx = 0; task_idx < num_cells; ++task_idx)
    cell_successor_offsets_[task_idx + 1] =
      cell_successor_offsets_[task_idx] +
      static_cast<std::uint32_t>(task_list[task_idx].successors.size());

  cell_successors_.resize(cell_successor_offsets_.back());
  for (std::size_t task_idx = 0; task_idx < num_cells; ++task_idx)
  {
    const auto& task = task_list[task_idx];
    std::copy(task.successors.begin(),
              task.successors.end(),
              cell_successors_.begin() + cell_successor_offsets_[task_idx]);
  }
}

void
CBCD_AngleSet::InitializeSweepState()
{
  const auto& initial_dependencies = cbc_spds_.GetInitialTaskDependencyCounts();
  std::copy(
    initial_dependencies.begin(), initial_dependencies.end(), remaining_cell_dependencies_.begin());
  batch_pipeline_.Reset();
  const auto& initial_ready_cells = cbc_spds_.GetInitialReadyTasks();
  auto& ready_cell_ids = batch_pipeline_.ReadyCells();
  ready_cell_ids.insert(
    ready_cell_ids.end(), initial_ready_cells.begin(), initial_ready_cells.end());
  num_completed_cells_ = 0;
  pending_reflecting_cells_ = following_angle_sets_.empty() ? 0 : num_reflecting_cells_;
}

CBCD_AngleSet::CompletedBatch
CBCD_AngleSet::RetireCompletedBatch()
{
  auto& completed_cell_ids = batch_pipeline_.cells[batch_pipeline_.launch_buffer];
  auto& ready_cell_ids = batch_pipeline_.ReadyCells();
  for (std::uint32_t i = 0; i < batch_pipeline_.launch_count; ++i)
  {
    const auto cell_local_id = completed_cell_ids[i];
    const auto succ_begin = cell_successor_offsets_[cell_local_id];
    const auto succ_end = cell_successor_offsets_[cell_local_id + 1];
    for (auto succ_i = succ_begin; succ_i < succ_end; ++succ_i)
    {
      if (--remaining_cell_dependencies_[cell_successors_[succ_i]] == 0)
        ready_cell_ids.push_back(cell_successors_[succ_i]);
    }

    if ((not following_angle_sets_.empty()) and (not followers_released_) and
        (has_outgoing_reflecting_face_[cell_local_id] != 0))
    {
      --pending_reflecting_cells_;
    }
  }

  num_completed_cells_ += batch_pipeline_.launch_count;
  const CompletedBatch completed{batch_pipeline_.launch_buffer, batch_pipeline_.launch_count};
  batch_pipeline_.launch_count = 0;
  return completed;
}

void
CBCD_AngleSet::AddReadyBatch(CBCDSweepChunk& sweep_chunk, const std::size_t worker_id)
{
  auto& ready_cell_ids = batch_pipeline_.ReadyCells();
  const auto launch_count = static_cast<std::uint32_t>(ready_cell_ids.size());
  batch_pipeline_.launch_buffer = batch_pipeline_.ready_buffer;
  batch_pipeline_.launch_count = launch_count;
  batch_pipeline_.ready_buffer = batch_pipeline_.AcquireFreeBuffer();
  batch_pipeline_.ReadyCells().clear();
  sweep_chunk.AddBatchSlot(worker_id, GetID(), ready_cell_ids.data(), launch_count);
}

void
CBCD_AngleSet::PublishCompletedBatch(const std::size_t worker_id, const CompletedBatch& completed)
{
  auto& completed_cell_ids = batch_pipeline_.cells[completed.buffer];
  cbcd_fluds_.PublishOutgoingPsi(
    *async_comm_, worker_id, GetID(), {completed_cell_ids.data(), std::size_t{completed.count}});
  completed_cell_ids.clear();
  batch_pipeline_.ReleaseBuffer(completed.buffer);
}

void
CBCD_AngleSet::TryReleaseFollowers()
{
  if (followers_released_)
    return;

  if (following_angle_sets_.empty())
  {
    followers_released_ = true;
    return;
  }

  if (pending_reflecting_cells_ != 0)
    return;

  // Publish reflecting-boundary writes before releasing follower dependencies.
  for (auto* angle_set : following_angle_sets_)
  {
    auto& cbcd_angle_set = dynamic_cast<CBCD_AngleSet&>(*angle_set);
    cbcd_angle_set.unresolved_sweep_dependencies_.fetch_sub(1, std::memory_order_release);
  }
  followers_released_ = true;
}

bool
CBCD_AngleSet::TryInitialize()
{
  if (sweep_initialized_)
    return false;
  if (unresolved_sweep_dependencies_.load(std::memory_order_acquire) != 0)
    return false;

  CALI_CXX_MARK_SCOPE("CBCD_AngleSet::TryInitialize");

  InitializeSweepState();
  sweep_initialized_ = true;
  return true;
}

void
CBCD_AngleSet::FinishCompletion(CBCDSweepChunk& sweep_chunk)
{
  CALI_CXX_MARK_SCOPE("CBCD_AngleSet::FinishCompletion");
  async_comm_->SignalAngleSetComplete(GetID());
  TryReleaseFollowers();
  executed_ = true;
  cbcd_fluds_.CopySaveAngularFluxFromDevice();
  cbcd_fluds_.CopySaveAngularFluxToDestinationPsi(
    sweep_chunk.GetProblem(), sweep_chunk.GetGroupset(), *this);
}

bool
CBCD_AngleSet::TryAdvanceOneStep(CBCDSweepChunk& sweep_chunk, const std::size_t worker_id)
{
  CALI_CXX_MARK_SCOPE("CBCD_AngleSet::TryAdvanceOneStep");

  if (executed_ or (not sweep_initialized_))
    return false;

  auto& ready_cell_ids = batch_pipeline_.ReadyCells();
  const bool kernel_completed =
    batch_pipeline_.HasKernelInFlight() and sweep_chunk.IsLaneLaunchCompleted(worker_id);
  const bool has_incoming = async_comm_->HasIncoming(GetID());

  // A sweep without local cells is drained before anything happens.
  if ((not kernel_completed) and ready_cell_ids.empty() and (not has_incoming) and
      (not IsSweepDrained()))
    return false;

  bool work_done = false;

  CompletedBatch completed;
  if (kernel_completed)
  {
    CALI_CXX_MARK_SCOPE("CBCD_AngleSet::RetireCompletedBatch");
    completed = RetireCompletedBatch();
    TryReleaseFollowers();
    work_done = true;
  }

  if (has_incoming)
  {
    CALI_CXX_MARK_SCOPE("CBCD_AngleSet::ProcessIncomingDependencies");
    work_done |= async_comm_->ProcessIncomingDependencies(
      GetID(),
      [&ready_cell_ids, this](std::uint32_t cell_local_id)
      {
        if (--remaining_cell_dependencies_[cell_local_id] == 0)
          ready_cell_ids.push_back(cell_local_id);
      });
  }

  if ((not batch_pipeline_.HasKernelInFlight()) and (not ready_cell_ids.empty()) and
      (not sweep_chunk.IsLaneInFlight(worker_id)))
  {
    CALI_CXX_MARK_SCOPE("CBCD_AngleSet::AddReadyBatch");
    AddReadyBatch(sweep_chunk, worker_id);
    work_done = true;
  }

  // Publish completed output after admitting ready work to the next worker launch.
  if (completed.count != 0)
  {
    CALI_CXX_MARK_SCOPE("CBCD_AngleSet::PublishCompletedBatch");
    PublishCompletedBatch(worker_id, completed);
    work_done = true;
  }

  if (IsSweepDrained())
  {
    FinishCompletion(sweep_chunk);
    return true;
  }

  return work_done;
}

AngleSetStatus
CBCD_AngleSet::AngleSetAdvance(SweepChunk& /*sweep_chunk*/, AngleSetStatus /*permission*/)
{
  OpenSnLogicalError("Device CBC angle sets are advanced only by the ASYNC_FIFO scheduler.");
}

void
CBCD_AngleSet::ResetSweepBuffers()
{
  batch_pipeline_.Reset();
  num_completed_cells_ = 0;
  pending_reflecting_cells_ = 0;
  sweep_initialized_ = false;
  followers_released_ = false;
  ResetSweepDependencies();
  executed_ = false;
}

} // namespace opensn
