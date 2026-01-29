// SPDX-FileCopyrightText: 2025 The OpenSn Authors <https://open-sn.github.io/opensn/>
// SPDX-License-Identifier: MIT

#include "modules/linear_boltzmann_solvers/discrete_ordinates_problem/sweep/scheduler/sweep_scheduler.h"
#include "modules/linear_boltzmann_solvers/discrete_ordinates_problem/sweep/spds/cbc.h"
#include "modules/linear_boltzmann_solvers/discrete_ordinates_problem/sweep/communicators/cbc_async_comm.h"
#include "modules/linear_boltzmann_solvers/discrete_ordinates_problem/sweep_chunks/cbcd_sweep_chunk.h"
#include "caribou/main.hpp"
#include "framework/logging/log.h"
#include "framework/runtime.h"
#include "caliper/cali.h"

namespace opensn
{

void
SweepScheduler::ScheduleAlgoAsyncFIFO(SweepChunk& sweep_chunk)
{
  CALI_CXX_MARK_SCOPE("SweepScheduler::ScheduleAlgoAsyncFIFO");

  CBCD_SweepChunk& cbcd_sweep_chunk = static_cast<CBCD_SweepChunk&>(sweep_chunk);
  const bool is_saving_angular_fluxes = cbcd_sweep_chunk.IsSavingAngularFluxes();

  cbcd_sweep_chunk.CopyPhiAndSrcToDevice();

  const auto& angle_sets = cbcd_sweep_chunk.GetAngleSets();
  const auto& fluds_list = cbcd_sweep_chunk.GetFLUDSList();
  const auto& streams_list = cbcd_sweep_chunk.GetStreamsList();

  const size_t num_angle_sets = angle_sets.size();
  std::vector<std::uint8_t> executed(num_angle_sets, 0);
  std::vector<std::uint8_t> boundary_data_set(num_angle_sets, 0);
  std::vector<std::uint8_t> kernel_in_flight(num_angle_sets, 0);
  std::vector<std::vector<Task*>> ready_tasks(num_angle_sets);
  std::vector<std::vector<Task*>> in_flight_tasks(num_angle_sets);

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
        fluds_list[i]->CopyOutgoingPsiBackToHost(cbcd_sweep_chunk, angle_sets[i]);
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
        kernel_in_flight[i] = 0;
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
      if (executed[i] or kernel_in_flight[i] or boundary_data_set[i])
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
        boundary_data_set[i] = 1;
        any_work_done = true;
      }
    }

    // Collect ready tasks and launch kernels (only if task dependencies changed)
    if (any_work_done)
    {
      for (size_t i = 0; i < num_angle_sets; ++i)
      {
        if (executed[i] or not boundary_data_set[i] or kernel_in_flight[i])
          continue;
        ready_tasks[i].clear();
        std::vector<std::uint64_t> ready_cell_ids;
        for (auto& task : angle_sets[i]->GetCurrentTaskList())
          if (task.num_dependencies == 0 and not task.completed)
          {
            ready_tasks[i].push_back(&task);
            ready_cell_ids.push_back(task.reference_id);
          }
        if (ready_tasks[i].empty())
          continue;
        fluds_list[i]->CopyIncomingNonlocalPsiToDevice(angle_sets[i], ready_cell_ids);
        cbcd_sweep_chunk.GPUSweep(*angle_sets[i], ready_cell_ids);
        fluds_list[i]->CopyOutgoingPsiFromDevice(angle_sets[i], ready_cell_ids);
        in_flight_tasks[i] = std::move(ready_tasks[i]);
        kernel_in_flight[i] = 1;
      }
    }

    // Check angleset completion
    for (size_t i = 0; i < num_angle_sets; ++i)
    {
      if (executed[i] or kernel_in_flight[i])
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
        executed[i] = 1;
        ++executed_anglesets;
        if (is_saving_angular_fluxes)
        {
          cbcd_sweep_chunk.CopySavedPsiFromDevice(*angle_sets[i]);
          streams_list[i]->add_callback(
            [&, i]() { cbcd_sweep_chunk.CopySavedPsiBackToHost(*angle_sets[i]); });
        }
      }
    }
  }

  cbcd_sweep_chunk.CopyOutflowAndPhiBackToHost();

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