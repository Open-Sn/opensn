// SPDX-FileCopyrightText: 2024 The OpenSn Authors <https://open-sn.github.io/opensn/>
// SPDX-License-Identifier: MIT

#include "modules/linear_boltzmann_solvers/discrete_ordinates_problem/sweep/angle_set/cbc_angle_set.h"
#include "modules/linear_boltzmann_solvers/discrete_ordinates_problem/sweep/communicators/cbc_async_comm.h"
#include "modules/linear_boltzmann_solvers/discrete_ordinates_problem/sweep/fluds/cbc_fluds.h"
#include "modules/linear_boltzmann_solvers/discrete_ordinates_problem/sweep/spds/cbc.h"
#include "modules/linear_boltzmann_solvers/discrete_ordinates_problem/sweep_chunks/sweep_chunk.h"
#include "framework/mesh/mesh_continuum/mesh_continuum.h"
#include "framework/data_types/range.h"
#include "framework/logging/log.h"
#include "framework/runtime.h"
#include "caliper/cali.h"
#include <cassert>

namespace opensn
{

CBC_AngleSet::CBC_AngleSet(size_t id,
                           unsigned int num_groups,
                           const SPDS& spds,
                           std::shared_ptr<FLUDS>& fluds,
                           const std::vector<size_t>& angle_indices,
                           std::map<uint64_t, std::shared_ptr<SweepBoundary>>& boundaries,
                           const MPICommunicatorSet& comm_set)
  : AngleSet(id, num_groups, spds, fluds, angle_indices, boundaries),
    cbc_spds_(dynamic_cast<const CBC_SPDS&>(spds_)),
    ready_tasks_(),
    async_comm_(id, *fluds, comm_set),
    cbc_fluds_(dynamic_cast<CBC_FLUDS&>(*fluds))
{

  const auto& task_list = cbc_spds_.GetTaskList();
  const auto num_tasks = task_list.size();
  initial_dependencies_.resize(num_tasks);
  remaining_dependencies_.resize(num_tasks);
  initial_ready_tasks_.reserve(num_tasks);
  ready_tasks_.reserve(num_tasks);

  for (std::uint32_t task_idx = 0; task_idx < num_tasks; ++task_idx)
  {
    const auto& task = task_list[task_idx];
    const auto num_dependencies = task.num_dependencies;
    initial_dependencies_[task_idx] = num_dependencies;
    if (num_dependencies == 0)
      initial_ready_tasks_.push_back(task_idx);
  }

  ResetTaskState();
}

AsynchronousCommunicator*
CBC_AngleSet::GetCommunicator()
{
  return static_cast<AsynchronousCommunicator*>(&async_comm_);
}

AngleSetStatus
CBC_AngleSet::AngleSetAdvance(SweepChunk& sweep_chunk, AngleSetStatus permission)
{
  CALI_CXX_MARK_SCOPE("CBC_AngleSet::AngleSetAdvance");

  if (executed_)
    return AngleSetStatus::FINISHED;

  const auto& task_list = cbc_spds_.GetTaskList();
  sweep_chunk.SetAngleSet(*this);

  const auto tasks_who_received_data = async_comm_.ReceiveData();

  for (const auto& task_number : tasks_who_received_data)
  {
    assert(remaining_dependencies_[task_number] > 0);
    if (--remaining_dependencies_[task_number] == 0)
      ready_tasks_.push_back(task_number);
  }

  async_comm_.SendData();

  // Check if boundaries allow for execution
  for (auto& [bid, boundary] : boundaries_)
    if (not boundary->CheckAnglesReadyStatus(angles_))
      return AngleSetStatus::NOT_FINISHED;

  while (not ready_tasks_.empty())
  {
    const auto task_idx = ready_tasks_.back();
    ready_tasks_.pop_back();
    const auto& cell_task = task_list[task_idx];

    sweep_chunk.SetCell(cell_task.cell_ptr, *this);
    sweep_chunk.Sweep(*this);

    for (const auto& local_task_num : cell_task.successors)
    {
      assert(remaining_dependencies_[local_task_num] > 0);
      if (--remaining_dependencies_[local_task_num] == 0)
        ready_tasks_.push_back(local_task_num);
    }

    ++num_completed_tasks_;
    async_comm_.SendData();
  }

  const bool all_tasks_completed = (num_completed_tasks_ == task_list.size());
  const bool all_messages_sent = async_comm_.SendData();

  if (all_tasks_completed and all_messages_sent)
  {
    // Update boundary readiness
    for (auto& [bid, boundary] : boundaries_)
      boundary->UpdateAnglesReadyStatus(angles_);
    executed_ = true;
    return AngleSetStatus::FINISHED;
  }

  return AngleSetStatus::NOT_FINISHED;
}

void
CBC_AngleSet::ResetSweepBuffers()
{
  ResetTaskState();
  async_comm_.Reset();
  fluds_->ClearLocalAndReceivePsi();
  executed_ = false;
}

void
CBC_AngleSet::ResetTaskState()
{
  std::copy(
    initial_dependencies_.begin(), initial_dependencies_.end(), remaining_dependencies_.begin());
  ready_tasks_ = initial_ready_tasks_;
  num_completed_tasks_ = 0;
}

const double*
CBC_AngleSet::PsiBoundary(uint64_t boundary_id,
                          unsigned int angle_num,
                          uint64_t cell_local_id,
                          unsigned int face_num,
                          unsigned int fi,
                          unsigned int g,
                          bool surface_source_active)
{
  if ((not boundaries_[boundary_id]->IsReflecting()) and (not surface_source_active))
    return boundaries_[boundary_id]->ZeroFlux(g);
  return boundaries_[boundary_id]->PsiIncoming(cell_local_id, face_num, fi, angle_num, g);
}

double*
CBC_AngleSet::PsiReflected(uint64_t boundary_id,
                           unsigned int angle_num,
                           uint64_t cell_local_id,
                           unsigned int face_num,
                           unsigned int fi)
{
  return boundaries_[boundary_id]->PsiOutgoing(cell_local_id, face_num, fi, angle_num);
}

} // namespace opensn
