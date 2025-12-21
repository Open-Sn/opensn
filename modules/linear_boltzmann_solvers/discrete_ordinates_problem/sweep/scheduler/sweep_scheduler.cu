// SPDX-FileCopyrightText: 2025 The OpenSn Authors <https://open-sn.github.io/opensn/>
// SPDX-License-Identifier: MIT

#include "modules/linear_boltzmann_solvers/discrete_ordinates_problem/sweep/scheduler/sweep_scheduler.h"
#include "modules/linear_boltzmann_solvers/discrete_ordinates_problem/sweep/spds/cbc.h"
#include "modules/linear_boltzmann_solvers/discrete_ordinates_problem/sweep/angle_set/cbc_angle_set.h"
#include "modules/linear_boltzmann_solvers/discrete_ordinates_problem/sweep/angle_set/cbc_angle_set_helpers.h"
#include "modules/linear_boltzmann_solvers/discrete_ordinates_problem/sweep/communicators/cbc_async_comm.h"
#include "modules/linear_boltzmann_solvers/discrete_ordinates_problem/sweep_chunks/cbc_sweep_chunk.h"
#include "modules/linear_boltzmann_solvers/discrete_ordinates_problem/sweep/fluds/cbcd_fluds.h"
#include "modules/linear_boltzmann_solvers/discrete_ordinates_problem/sweep/boundary/reflecting_boundary.h"
#include "caribou/caribou.h"
#include "caliper/cali.h"

namespace opensn
{
void
SweepScheduler::DeviceScheduleAlgoFIFO(SweepChunk& sweep_chunk)
{
  CALI_CXX_MARK_SCOPE("SweepScheduler::DeviceScheduleAlgoFIFO");

  CBCSweepChunk& cbc_sweep_chunk = dynamic_cast<CBCSweepChunk&>(sweep_chunk);

  // Copy phi and source to device
  cbc_sweep_chunk.CopyPhiAndSrcToDevice();

  std::vector<CBC_AngleSet*> angle_sets;
  for (auto& angle_set_group : angle_agg_.angle_set_groups)
    for (auto& angle_set : angle_set_group.GetAngleSets())
      angle_sets.push_back(dynamic_cast<CBC_AngleSet*>(angle_set.get()));

  const size_t num_angle_sets = angle_sets.size();

  std::vector<bool> executed(num_angle_sets, false);
  std::vector<bool> boundary_data_set(num_angle_sets, false);
  std::vector<bool> kernel_in_flight(num_angle_sets, false);
  std::vector<std::vector<Task*>> ready_tasks(num_angle_sets);
  std::vector<std::vector<Task*>> in_flight_tasks(num_angle_sets);

  for (auto* angle_set : angle_sets)
  {
    auto& current_task_list = angle_set->GetCurrentTaskList();
    if (current_task_list.empty())
      current_task_list = dynamic_cast<const CBC_SPDS&>(angle_set->GetSPDS()).GetTaskList();
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

      crb::Stream& stream = GetCBCAngleSetStream(*angle_sets[i]);

      // Check if the kernel is done
      if (stream.is_completed())
      {
        // Copy back outgoing boundary (reflecting) and non-local psi
        auto& fluds = dynamic_cast<CBCD_FLUDS&>(angle_sets[i]->GetFLUDS());
        fluds.CopyOutgoingBoundaryPsiToHost(sweep_chunk, in_flight_tasks[i]);
        fluds.CopyOutgoingNonLocalPsiToHost(sweep_chunk, in_flight_tasks[i]);

        // Update task dependencies
        auto& current_task_list = angle_sets[i]->GetCurrentTaskList();
        for (auto* task : in_flight_tasks[i])
        {
          for (uint64_t succ : task->successors)
            --current_task_list[succ].num_dependencies;
          task->completed = true;
        }

        // Send MPI data
        auto* comm = dynamic_cast<CBC_ASynchronousCommunicator*>(angle_sets[i]->GetCommunicator());
        comm->SendData();

        in_flight_tasks[i].clear();
        kernel_in_flight[i] = false;
        any_work_done = true;
      }
    }

    // Receive and send MPI data
    for (size_t i = 0; i < num_angle_sets; ++i)
    {
      if (executed[i])
        continue;

      auto* comm = dynamic_cast<CBC_ASynchronousCommunicator*>(angle_sets[i]->GetCommunicator());
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
        dynamic_cast<CBCD_FLUDS&>(as->GetFLUDS()).CopyIncomingBoundaryPsiToDevice(sweep_chunk);
        boundary_data_set[i] = true;
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
        for (auto& task : angle_sets[i]->GetCurrentTaskList())
          if (task.num_dependencies == 0 and not task.completed)
            ready_tasks[i].push_back(&task);

        if (ready_tasks[i].empty())
          continue;

        auto& fluds = dynamic_cast<CBCD_FLUDS&>(angle_sets[i]->GetFLUDS());
        fluds.CopyIncomingNonLocalPsiToDevice(sweep_chunk, ready_tasks[i]);

        cbc_sweep_chunk.GPUSweep(*angle_sets[i], ready_tasks[i]);

        in_flight_tasks[i] = std::move(ready_tasks[i]);
        kernel_in_flight[i] = true;
      }
    }

    // Check angleset completion
    for (size_t i = 0; i < num_angle_sets; ++i)
    {
      if (executed[i] or kernel_in_flight[i])
        continue;

      auto& current_task_list = angle_sets[i]->GetCurrentTaskList();
      auto* comm = dynamic_cast<CBC_ASynchronousCommunicator*>(angle_sets[i]->GetCommunicator());

      bool all_done = std::all_of(current_task_list.begin(),
                                  current_task_list.end(),
                                  [](const Task& t) { return t.completed; });

      if (all_done and comm->SendData())
      {
        for (auto& [bid, boundary] : angle_sets[i]->GetBoundaries())
          boundary->UpdateAnglesReadyStatus(angle_sets[i]->GetAngleIndices());
        executed[i] = true;
        ++executed_anglesets;
      }
    }
  }

  // Copy back phi and outflow from device
  cbc_sweep_chunk.CopyOutflowAndPhiFromDevice();

  opensn::mpi_comm.barrier();

  bool received_delayed_data = false;
  while (not received_delayed_data)
  {
    received_delayed_data = true;

    for (auto& angle_set_group : angle_agg_.angle_set_groups)
      for (auto& angle_set : angle_set_group.GetAngleSets())
      {
        if (angle_set->FlushSendBuffers() == AngleSetStatus::MESSAGES_PENDING)
          received_delayed_data = false;

        if (not angle_set->ReceiveDelayedData())
          received_delayed_data = false;
      }
  }

  for (auto& angle_set_group : angle_agg_.angle_set_groups)
    for (auto& angle_set : angle_set_group.GetAngleSets())
      angle_set->ResetSweepBuffers();

  for (const auto& [bid, bndry] : angle_agg_.GetSimBoundaries())
  {
    if (bndry->GetType() == LBSBoundaryType::REFLECTING)
    {
      auto rbndry = std::static_pointer_cast<ReflectingBoundary>(bndry);
      rbndry->ResetAnglesReadyStatus();
    }
  }
}
} // namespace opensn