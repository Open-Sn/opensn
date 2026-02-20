// SPDX-FileCopyrightText: 2026 The OpenSn Authors <https://open-sn.github.io/opensn/>
// SPDX-License-Identifier: MIT

#include "modules/linear_boltzmann_solvers/discrete_ordinates_problem/sweep/scheduler/sweep_scheduler.h"
#include "modules/linear_boltzmann_solvers/discrete_ordinates_problem/sweep/angle_set/aahd_angle_set.h"
#include "modules/linear_boltzmann_solvers/discrete_ordinates_problem/sweep_chunks/aahd_sweep_chunk.h"
#include "modules/linear_boltzmann_solvers/discrete_ordinates_problem/sweep/spds/cbc.h"
#include "modules/linear_boltzmann_solvers/discrete_ordinates_problem/sweep/communicators/cbc_async_comm.h"
#include "modules/linear_boltzmann_solvers/discrete_ordinates_problem/sweep_chunks/cbcd_sweep_chunk.h"
#include "modules/linear_boltzmann_solvers/discrete_ordinates_problem/discrete_ordinates_problem.h"
#include "caribou/main.hpp"
#include "caliper/cali.h"
#include <thread>
#include <vector>

namespace opensn
{

void
SweepScheduler::ScheduleAlgoAAO(SweepChunk& sweep_chunk)
{
  CALI_CXX_MARK_SCOPE("SweepScheduler::ScheduleAlgoAAO");

  // copy phi and src moments to device
  auto aah_sweep_chunk = static_cast<AAHDSweepChunk&>(sweep_chunk);
  aah_sweep_chunk.GetProblem().CopyPhiAndSrcToDevice();

  // allocate threads for each angleset
  std::size_t num_anglesets = angle_agg_.GetNumAngleSets();

  // set the latches for all anglesets
  for (auto& angle_set : angle_agg_)
  {
    auto aahd_angle_set = static_cast<AAHD_AngleSet*>(angle_set.get());
    aahd_angle_set->SetStartingLatch();
  }

  // launch threads
  pool_.run(
    [this, &sweep_chunk](std::size_t i)
    {
      auto* aahd = static_cast<AAHD_AngleSet*>(angle_agg_[i].get());
      aahd->AngleSetAdvance(sweep_chunk, AngleSetStatus::EXECUTE);
    });

  // copy phi and outflow data back to host
  aah_sweep_chunk.GetProblem().CopyPhiAndOutflowBackToHost();

  // reset all anglesets
  for (auto& angle_set : angle_agg_)
    angle_set->ResetSweepBuffers();
}

void
SweepScheduler::ScheduleAlgoAsyncFIFO(SweepChunk& sweep_chunk)
{
  CALI_CXX_MARK_SCOPE("SweepScheduler::ScheduleAlgoAsyncFIFO");

  auto& cbcd_sweep_chunk = static_cast<CBCDSweepChunk&>(sweep_chunk);
  // Copy phi and source moments to device
  cbcd_sweep_chunk.GetProblem().CopyPhiAndSrcToDevice();

  const auto& angle_sets = cbcd_sweep_chunk.GetAngleSets();
  const auto& fluds_list = cbcd_sweep_chunk.GetFLUDSList();
  const auto& streams_list = cbcd_sweep_chunk.GetStreamsList();

  const size_t num_angle_sets = angle_sets.size();
  std::vector<bool> executed(num_angle_sets, 0);
  std::vector<bool> boundary_data_set(num_angle_sets, 0);
  std::vector<bool> kernel_in_flight(num_angle_sets, 0);
  std::vector<std::vector<Task*>> ready_tasks(num_angle_sets);
  std::vector<std::vector<std::uint64_t>> ready_cell_ids(num_angle_sets);
  std::vector<std::vector<Task*>> in_flight_tasks(num_angle_sets);
  std::vector<std::vector<std::uint64_t>> in_flight_cell_ids(num_angle_sets);

  for (auto* angle_set : angle_sets)
  {
    auto& current_task_list = angle_set->GetCurrentTaskList();
    if (current_task_list.empty())
      current_task_list = static_cast<const CBC_SPDS&>(angle_set->GetSPDS()).GetTaskList();
  }

  size_t executed_anglesets = 0;
  while (executed_anglesets < num_angle_sets)
  {
    bool any_work_done = false;

    // Poll completed kernels
    for (size_t i = 0; i < num_angle_sets; ++i)
    {
      if (not kernel_in_flight[i])
        continue;
      // Check if the kernel is done
      if (streams_list[i]->is_completed())
      {
        // Copy back outgoing (reflecting) boundary and non-local psi
        fluds_list[i]->CopyOutgoingPsiBackToHost(
          cbcd_sweep_chunk, angle_sets[i], in_flight_cell_ids[i]);
        // Update task dependencies
        auto& current_task_list = angle_sets[i]->GetCurrentTaskList();
        for (auto* task : in_flight_tasks[i])
        {
          for (uint64_t succ : task->successors)
            --current_task_list[succ].num_dependencies;
          task->completed = true;
        }
        // Send MPI data
        auto* comm = static_cast<CBC_AsynchronousCommunicator*>(angle_sets[i]->GetCommunicator());
        comm->SendData();
        in_flight_tasks[i].clear();
        in_flight_cell_ids[i].clear();
        kernel_in_flight[i] = false;
        any_work_done = true;
      }
    }

    // Receive and send MPI data
    for (size_t i = 0; i < num_angle_sets; ++i)
    {
      if (executed[i])
        continue;
      auto* comm = static_cast<CBC_AsynchronousCommunicator*>(angle_sets[i]->GetCommunicator());
      auto& current_task_list = angle_sets[i]->GetCurrentTaskList();
      auto received = comm->ReceiveData();
      if (not received.empty())
      {
        for (uint64_t t : received)
          --current_task_list[t].num_dependencies;
        any_work_done = true;
      }
      comm->SendData();
    }

    // Set boundary data
    for (size_t i = 0; i < num_angle_sets; ++i)
    {
      if (executed[i] or boundary_data_set[i] or kernel_in_flight[i])
        continue;
      auto* as = angle_sets[i];
      bool boundaries_ready = true;
      for (auto& [bid, boundary] : as->GetBoundaries())
      {
        if (not boundary->CheckAnglesReadyStatus(as->GetAngleIndices()))
        {
          boundaries_ready = false;
          break;
        }
      }
      if (boundaries_ready)
      {
        fluds_list[i]->CopyIncomingBoundaryPsiToDevice(cbcd_sweep_chunk, angle_sets[i]);
        boundary_data_set[i] = true;
        any_work_done = true;
      }
    }

    // Collect ready tasks and launch kernels (only if task dependencies changed)
    if (any_work_done)
    {
      for (size_t i = 0; i < num_angle_sets; ++i)
      {
        if (executed[i] or (not boundary_data_set[i]) or kernel_in_flight[i])
          continue;
        ready_tasks[i].clear();
        ready_cell_ids[i].clear();
        for (auto& task : angle_sets[i]->GetCurrentTaskList())
          if (task.num_dependencies == 0 and not task.completed)
          {
            ready_tasks[i].push_back(&task);
            ready_cell_ids[i].push_back(task.reference_id);
          }
        if (ready_tasks[i].empty())
          continue;
        fluds_list[i]->CopyIncomingNonlocalPsiToDevice(angle_sets[i], ready_cell_ids[i]);
        cbcd_sweep_chunk.GPUSweep(*angle_sets[i], ready_cell_ids[i]);
        in_flight_tasks[i] = std::move(ready_tasks[i]);
        in_flight_cell_ids[i] = std::move(ready_cell_ids[i]);
        kernel_in_flight[i] = true;
      }
    }

    // Check angleset completion
    for (size_t i = 0; i < num_angle_sets; ++i)
    {
      if (executed[i] or (not boundary_data_set[i]) or kernel_in_flight[i])
        continue;
      auto& current_task_list = angle_sets[i]->GetCurrentTaskList();
      auto* comm = static_cast<CBC_AsynchronousCommunicator*>(angle_sets[i]->GetCommunicator());
      bool all_done = std::all_of(current_task_list.begin(),
                                  current_task_list.end(),
                                  [](const Task& t) { return t.completed; });
      if (all_done and comm->SendData())
      {
        for (auto& [bid, boundary] : angle_sets[i]->GetBoundaries())
          boundary->UpdateAnglesReadyStatus(angle_sets[i]->GetAngleIndices());
        executed[i] = true;
        ++executed_anglesets;
        fluds_list[i]->CopySavedPsiFromDevice();
        streams_list[i]->add_callback(
          [&, i]()
          { fluds_list[i]->CopySavedPsiToDestinationPsi(cbcd_sweep_chunk, angle_sets[i]); });
      }
    }
  }

  /// Copy phi and outflow data back to host
  cbcd_sweep_chunk.GetProblem().CopyPhiAndOutflowBackToHost();

  // Receive delayed data
  opensn::mpi_comm.barrier();
  bool received_delayed_data = false;
  while (not received_delayed_data)
  {
    received_delayed_data = true;

    for (auto& angle_set : angle_sets)
    {
      if (angle_set->FlushSendBuffers() == AngleSetStatus::MESSAGES_PENDING)
        received_delayed_data = false;

      if (not angle_set->ReceiveDelayedData())
        received_delayed_data = false;
    }
  }

  // Reset all
  for (auto& angle_set : angle_sets)
    angle_set->ResetSweepBuffers();

  for (const auto& [bid, bndry] : angle_agg_.GetSimBoundaries())
    bndry->ResetAnglesReadyStatus();
}

} // namespace opensn
