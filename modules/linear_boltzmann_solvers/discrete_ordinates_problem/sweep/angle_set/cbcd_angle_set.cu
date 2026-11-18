// SPDX-FileCopyrightText: 2026 The OpenSn Authors <https://open-sn.github.io/opensn/>
// SPDX-License-Identifier: MIT

#include "modules/linear_boltzmann_solvers/discrete_ordinates_problem/sweep/angle_set/cbcd_angle_set.h"
#include "modules/linear_boltzmann_solvers/discrete_ordinates_problem/sweep/angle_aggregation/angle_aggregation.h"
#include "modules/linear_boltzmann_solvers/discrete_ordinates_problem/sweep/fluds/cbcd_fluds.h"
#include "modules/linear_boltzmann_solvers/discrete_ordinates_problem/sweep/spds/cbc.h"
#include "modules/linear_boltzmann_solvers/discrete_ordinates_problem/sweep_chunks/cbcd_sweep_chunk.h"
#include "modules/linear_boltzmann_solvers/discrete_ordinates_problem/sweep/boundary/sweep_boundary.h"
#include "framework/utils/error.h"
#include "caliper/cali.h"
#include <algorithm>
#include <bit>
#include <cassert>
#include <limits>
#include <map>
#include <set>
#include <stdexcept>
#include <string>
#include <utility>

namespace opensn
{

CBCD_AngleSet::CBCD_AngleSet(std::size_t id,
                             const LBSGroupset& groupset,
                             const SPDS& spds,
                             std::shared_ptr<FLUDS>& fluds,
                             const std::vector<std::size_t>& angle_indices,
                             std::map<std::uint64_t, std::shared_ptr<SweepBoundary>>& boundaries,
                             const SweepCommunicator& sweep_communicator)
  : AngleSet(id, groupset, spds, fluds, angle_indices, boundaries),
    cbc_spds_(dynamic_cast<const CBC_SPDS&>(spds)),
    sweep_communicator_(sweep_communicator),
    cbcd_fluds_(dynamic_cast<CBCD_FLUDS&>(*fluds_)),
    device_angle_indices_(angles_.size())
{
  const auto num_cells = cbc_spds_.GetTaskList().size();
  ready_cells_.reserve(num_cells);
  batch_cells_.reserve(num_cells);
  BuildCellTaskGraph();
}

void
CBCD_AngleSet::SyncDeviceAngleIndices()
{
  crb::copy(device_angle_indices_, angles_, angles_.size());
}

void
CBCD_AngleSet::BuildReflectingDependencies(const std::vector<CBCD_AngleSet*>& angle_sets,
                                           const AngleAggregation& angle_agg)
{
  for (auto* angle_set : angle_sets)
  {
    angle_set->follower_queues_.clear();
    angle_set->cell_followers_.clear();
    angle_set->leader_queues_.clear();
    angle_set->leader_dependent_cells_.clear();
  }

  const auto ByID = [](const AngleSet* lhs, const AngleSet* rhs)
  { return lhs->GetID() < rhs->GetID(); };
  for (auto* leader : angle_sets)
  {
    if (leader->following_angle_sets_.empty())
      continue;

    std::map<std::uint64_t, std::set<AngleSet*>> followers_by_boundary;
    std::set<AngleSet*> reflecting_followers;
    for (const auto& [boundary_id, boundary] : leader->boundaries_)
    {
      if (not boundary->IsReflecting())
        continue;
      std::set<AngleSet*> followers;
      boundary->GetFollowingAngleSets(leader->GetGroupsetID(), followers, angle_agg, *leader);
      if (followers.empty())
        continue;
      reflecting_followers.insert(followers.begin(), followers.end());
      followers_by_boundary.emplace(boundary_id, std::move(followers));
    }
    OpenSnLogicalErrorIf(
      reflecting_followers != std::set<AngleSet*>(leader->following_angle_sets_.begin(),
                                                  leader->following_angle_sets_.end()) or
        reflecting_followers.contains(leader),
      "Device CBC: the followers of angle set " + std::to_string(leader->GetID()) +
        " are not those of its reflecting boundaries.");

    std::vector<CBCD_AngleSet*> followers;
    followers.reserve(reflecting_followers.size());
    for (auto* follower : reflecting_followers)
      followers.push_back(&dynamic_cast<CBCD_AngleSet&>(*follower));
    std::sort(followers.begin(), followers.end(), ByID);
    if (followers.size() > std::numeric_limits<FollowerMask>::digits)
      throw std::runtime_error(
        "Device CBC: angle set " + std::to_string(leader->GetID()) + " has " +
        std::to_string(followers.size()) + " reflecting follower angle sets; at most " +
        std::to_string(std::numeric_limits<FollowerMask>::digits) + " are supported.");
    const auto& common_data = leader->cbcd_fluds_.GetCommonData();
    auto& cell_followers = leader->cell_followers_;
    for (const auto& face : common_data.GetOutgoingBoundaryFaces())
      if (const auto it = followers_by_boundary.find(face.boundary_id);
          it != followers_by_boundary.end())
      {
        if (cell_followers.empty())
          cell_followers.assign(common_data.GetNumLocalCells(), 0);
        auto& mask = cell_followers[face.cell_local_id];
        for (auto* follower : it->second)
        {
          const auto bit = std::lower_bound(followers.begin(), followers.end(), follower, ByID) -
                           followers.begin();
          mask |= FollowerMask{1} << bit;
        }
      }
    if (cell_followers.empty())
      continue;

    std::vector<std::size_t> num_follower_cells(followers.size(), 0);
    for (std::uint32_t cell_local_id = 0; cell_local_id < cell_followers.size(); ++cell_local_id)
      for (auto mask = cell_followers[cell_local_id]; mask != 0; mask &= mask - 1)
      {
        const auto i = static_cast<std::size_t>(std::countr_zero(mask));
        ++num_follower_cells[i];
        followers[i]->leader_dependent_cells_.push_back(cell_local_id);
      }
    for (std::size_t i = 0; i < followers.size(); ++i)
    {
      leader->follower_queues_.push_back(
        std::make_unique<LockFreeSPSCQueue<std::uint32_t>>(num_follower_cells[i]));
      if (num_follower_cells[i] != 0)
        followers[i]->leader_queues_.push_back(leader->follower_queues_.back().get());
    }
  }
}

void
CBCD_AngleSet::BuildCellTaskGraph()
{
  cell_successor_offsets_.assign(1, 0);
  cell_successors_.clear();
  for (const auto& task : cbc_spds_.GetTaskList())
  {
    cell_successors_.insert(cell_successors_.end(), task.successors.begin(), task.successors.end());
    cell_successor_offsets_.push_back(static_cast<std::uint32_t>(cell_successors_.size()));
  }
}

void
CBCD_AngleSet::ResetSweepBuffers()
{
  for (auto& queue : follower_queues_)
    queue->Reset();
  sweep_initialized_ = false;
  executed_ = false;
}

void
CBCD_AngleSet::InitializeSweepState()
{
  remaining_cell_dependencies_ = cbc_spds_.GetInitialTaskDependencyCounts();
  for (const auto cell_local_id : leader_dependent_cells_)
    ++remaining_cell_dependencies_[cell_local_id];
  ready_cells_.clear();
  for (const auto cell_local_id : cbc_spds_.GetInitialReadyTasks())
    if (remaining_cell_dependencies_[cell_local_id] == 0)
      ready_cells_.push_back(cell_local_id);
  batch_cells_.clear();
  batch_successors_released_ = false;
  num_completed_cells_ = 0;
}

void
CBCD_AngleSet::ResolveCellDependencies(std::span<const std::uint32_t> cell_local_ids)
{
  for (const auto cell_local_id : cell_local_ids)
    if (--remaining_cell_dependencies_[cell_local_id] == 0)
      ready_cells_.push_back(cell_local_id);
}

void
CBCD_AngleSet::ReleaseBatchSuccessors()
{
  for (const auto cell_local_id : batch_cells_)
  {
    const auto succ_begin = cell_successor_offsets_[cell_local_id];
    const auto succ_end = cell_successor_offsets_[cell_local_id + 1];
    for (auto succ_i = succ_begin; succ_i < succ_end; ++succ_i)
      if (--remaining_cell_dependencies_[cell_successors_[succ_i]] == 0)
        ready_cells_.push_back(cell_successors_[succ_i]);
  }
}

void
CBCD_AngleSet::CompleteBatch()
{
  if (not cell_followers_.empty())
    for (const auto cell_local_id : batch_cells_)
      for (auto mask = cell_followers_[cell_local_id]; mask != 0; mask &= mask - 1)
        follower_queues_[std::countr_zero(mask)]->Push(cell_local_id);
  num_completed_cells_ += batch_cells_.size();
}

void
CBCD_AngleSet::FinishCompletion(CBCDSweepChunk& sweep_chunk)
{
  CALI_CXX_MARK_SCOPE("CBCD_AngleSet::FinishCompletion");
  executed_ = true;
  cbcd_fluds_.CopySaveAngularFluxToDestinationPsi(
    sweep_chunk.GetProblem(), sweep_chunk.GetGroupset(), *this);
}

bool
CBCD_AngleSet::TryAdvanceOneStep(CBCDSweepChunk& sweep_chunk, const std::size_t worker_id)
{
  assert(not executed_);
  bool work_done = false;
  if (not sweep_initialized_)
  {
    CALI_CXX_MARK_SCOPE("CBCD_AngleSet::InitializeSweepState");
    InitializeSweepState();
    sweep_initialized_ = true;
    work_done = true;
  }

  // The visit that adds a batch submits it before this angle set advances again.
  if (not batch_cells_.empty() and not batch_successors_released_)
  {
    ReleaseBatchSuccessors();
    batch_successors_released_ = true;
    work_done = true;
  }

  if (not batch_cells_.empty() and sweep_chunk.IsLaneLaunchCompleted(worker_id))
  {
    CompleteBatch();
    batch_cells_.clear();
    work_done = true;
  }

  if (const auto cells = async_comm_->TakeIncomingDependencies(GetID()); not cells.empty())
  {
    ResolveCellDependencies(cells);
    work_done = true;
  }

  for (auto* queue : leader_queues_)
    if (const auto cells = queue->Take(); not cells.empty())
    {
      ResolveCellDependencies(cells);
      work_done = true;
    }

  if (batch_cells_.empty() and (not ready_cells_.empty()) and
      (not sweep_chunk.IsLaneInFlight(worker_id)))
  {
    std::swap(ready_cells_, batch_cells_);
    batch_successors_released_ = false;
    sweep_chunk.AddBatchSlot(
      worker_id, GetID(), batch_cells_.data(), static_cast<std::uint32_t>(batch_cells_.size()));
    work_done = true;
  }

  // Every pending batch contains at least one cell that has yet to complete.
  if (num_completed_cells_ == remaining_cell_dependencies_.size())
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

} // namespace opensn
