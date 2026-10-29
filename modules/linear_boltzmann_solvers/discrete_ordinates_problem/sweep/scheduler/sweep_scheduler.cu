// SPDX-FileCopyrightText: 2026 The OpenSn Authors <https://open-sn.github.io/opensn/>
// SPDX-License-Identifier: MIT

#include "modules/linear_boltzmann_solvers/discrete_ordinates_problem/sweep/scheduler/sweep_scheduler.h"
#include "modules/linear_boltzmann_solvers/discrete_ordinates_problem/sweep/angle_set/aahd_angle_set.h"
#include "modules/linear_boltzmann_solvers/discrete_ordinates_problem/sweep_chunks/aahd_sweep_chunk.h"
#include "modules/linear_boltzmann_solvers/discrete_ordinates_problem/sweep/spds/cbc.h"
#include "modules/linear_boltzmann_solvers/discrete_ordinates_problem/sweep/communicators/cbcd_async_comm.h"
#include "modules/linear_boltzmann_solvers/discrete_ordinates_problem/sweep_chunks/cbcd_sweep_chunk.h"
#include "modules/linear_boltzmann_solvers/discrete_ordinates_problem/discrete_ordinates_problem.h"
#include "caribou/main.hpp"
#include <algorithm>
#include <thread>
#include <vector>

namespace opensn
{

void
SweepScheduler::ScheduleAlgoAAO(SweepChunk& sweep_chunk)
{
  // copy phi and src moments to device
  auto& aah_sweep_chunk = dynamic_cast<AAHDSweepChunk&>(sweep_chunk);
  int groupset_id = angle_agg_.GetGroupsetID();
  aah_sweep_chunk.GetProblem().CopyPhiAndSrcToDevice();
  aah_sweep_chunk.GetProblem().TransferDeviceBoundaryData(groupset_id, true);

  // reset dependency counters and pre-post receives for all anglesets
  std::size_t num_anglesets = angle_agg_.GetNumAngleSets();
  execution_order_.resize(num_anglesets);
  for (std::size_t i = 0; i < num_anglesets; ++i)
  {
    auto* aahd_angle_set = dynamic_cast<AAHD_AngleSet*>(angle_agg_[i].get());
    aahd_angle_set->ResetDependencyCounter();
    aahd_angle_set->PrepostReceives();
    execution_order_[i] = i;
  }

  // assign sweep task to thread pool (but not execution yet)
  pool_.AssignTask(
    [this, &sweep_chunk](std::size_t i)
    {
      auto* aahd = dynamic_cast<AAHD_AngleSet*>(angle_agg_[i].get());
      aahd->AngleSetAdvance(sweep_chunk, AngleSetStatus::EXECUTE);
    });

  // poll for readiness and launch threads
  while (!execution_order_.empty())
  {
    for (auto it = execution_order_.begin(); it != execution_order_.end();)
    {
      auto* angle_set = dynamic_cast<AAHD_AngleSet*>(angle_agg_[*it].get());
      if (angle_set->IsReady())
      {
        pool_.Run(*it);
        std::swap(*it, execution_order_.back());
        execution_order_.pop_back();
      }
      else
      {
        ++it;
      }
    }
  }
  pool_.WaitAll();

  // wait for sends and receive of delayed data
  for (auto& angle_set : angle_agg_)
  {
    auto* aahd_angle_set = dynamic_cast<AAHD_AngleSet*>(angle_set.get());
    aahd_angle_set->WaitForDownstreamAndDelayed();
  }

  // copy boundary, phi and outflow data back to host
  aah_sweep_chunk.GetProblem().TransferDeviceBoundaryData(groupset_id, false);
  aah_sweep_chunk.GetProblem().CopyPhiAndOutflowBackToHost();

  // reset all anglesets
  for (auto& angle_set : angle_agg_)
    angle_set->ResetSweepBuffers();
}

void
SweepScheduler::ScheduleAlgoAAOFused(SweepChunk& sweep_chunk)
{
  auto& aah_sweep_chunk = dynamic_cast<AAHDSweepChunk&>(sweep_chunk);
  int groupset_id = angle_agg_.GetGroupsetID();
  aah_sweep_chunk.GetProblem().CopyPhiAndSrcToDevice();
  aah_sweep_chunk.GetProblem().TransferDeviceBoundaryData(groupset_id, true);

  // Every lane's stream is idle here. The delayed-psi uploads in PrepostReceives go on lane 0. A
  // run may use another lane, so lane 0 is synchronized before any run starts.
  std::size_t num_anglesets = angle_agg_.GetNumAngleSets();
  std::vector<AAHD_AngleSet*> all_anglesets(num_anglesets);
  crb::Stream& setup_stream = aah_sweep_chunk.GetLaneStream(0);
  for (std::size_t i = 0; i < num_anglesets; ++i)
  {
    auto* aahd_angle_set = dynamic_cast<AAHD_AngleSet*>(angle_agg_[i].get());
    aahd_angle_set->ResetDependencyCounter();
    aahd_angle_set->SetFLUDSStream(setup_stream);
    aahd_angle_set->PrepostReceives();
    all_anglesets[i] = aahd_angle_set;
  }
  setup_stream.synchronize();

  // Maximum angle sets per fused run. Each member is one grid-z slice, so this bounds the launch
  // well below the CUDA/HIP grid-z limit (65535) and bounds the per-run slot array at
  // max_batch_width * levels entries.
  constexpr std::size_t max_batch_width = 256;
  const std::size_t num_lanes = aah_sweep_chunk.NumFusedLanes();
  std::vector<std::vector<AAHD_AngleSet*>> lane_runs(num_lanes);
  // Angle sets that have not started, in angle-set order. Only these are polled for readiness.
  std::vector<AAHD_AngleSet*> pending = all_anglesets;
  std::size_t num_finished = 0;

  // Runs are pipelined. Up to num_lanes runs are in flight, each on its own stream. When a lane is
  // idle, its run's members are completed. This releases their followers and sends their
  // downstream angular fluxes. The angle sets that are ready then start on a free lane. They do
  // not wait for the other runs. A run's members are fixed when it starts, so the run needs no
  // host interaction until it is done.
  while (num_finished < num_anglesets)
  {
    for (std::size_t lane = 0; lane < num_lanes; ++lane)
    {
      auto& run = lane_runs[lane];
      if (run.empty() or not aah_sweep_chunk.IsLaneIdle(lane))
        continue;
      for (auto* as : run)
      {
        as->FinishCompletion(sweep_chunk);
        ++num_finished;
      }
      run.clear();
    }

    for (std::size_t lane = 0; lane < num_lanes; ++lane)
    {
      auto& run = lane_runs[lane];
      if (not run.empty())
        continue;
      crb::Stream& stream = aah_sweep_chunk.GetLaneStream(lane);
      auto still_pending = pending.begin();
      for (auto* as : pending)
      {
        if (run.size() < max_batch_width and as->IsReady())
        {
          as->PrepareForBatch(sweep_chunk, stream);
          run.push_back(as);
        }
        else
          *still_pending++ = as;
      }
      pending.erase(still_pending, pending.end());
      if (run.empty())
        break; // nothing else is ready yet
      aah_sweep_chunk.SubmitBatchRun(lane, run);
    }
  }

  // Wait for sends and receive of delayed data
  for (auto& angle_set : angle_agg_)
  {
    auto* aahd_angle_set = dynamic_cast<AAHD_AngleSet*>(angle_set.get());
    aahd_angle_set->WaitForDownstreamAndDelayed();
  }

  // Copy boundary, phi and outflow data back to host
  aah_sweep_chunk.GetProblem().TransferDeviceBoundaryData(groupset_id, false);
  aah_sweep_chunk.GetProblem().CopyPhiAndOutflowBackToHost();

  // Reset all angle sets
  for (auto& angle_set : angle_agg_)
    angle_set->ResetSweepBuffers();
}

void
SweepScheduler::ScheduleAlgoAsyncFIFO(SweepChunk& sweep_chunk)
{
  // Reset dependency counter
  for (auto& angle_set : angle_agg_)
    angle_set->ResetDependencyCounter();

  auto& cbcd_sweep_chunk = dynamic_cast<CBCDSweepChunk&>(sweep_chunk);
  // Copy phi and source moments to device
  cbcd_sweep_chunk.GetProblem().CopyPhiAndSrcToDevice();

  const auto& angle_sets = cbcd_sweep_chunk.GetAngleSets();
  const auto& fluds_list = cbcd_sweep_chunk.GetFLUDS();
  auto& streams_list = cbcd_sweep_chunk.GetStreams();

  const size_t num_angle_sets = angle_sets.size();
  std::vector<bool> executed(num_angle_sets, false);
  std::vector<bool> boundary_data_set(num_angle_sets, false);
  std::vector<bool> kernel_in_flight(num_angle_sets, false);
  std::vector<std::vector<Task*>> ready_queues(num_angle_sets);
  std::vector<size_t> num_completed_tasks(num_angle_sets, 0);
  std::vector<std::vector<Task*>> ready_tasks(num_angle_sets);
  std::vector<std::vector<std::uint32_t>> ready_cell_ids(num_angle_sets);
  std::vector<std::vector<Task*>> in_flight_tasks(num_angle_sets);
  std::vector<std::vector<std::uint32_t>> in_flight_cell_ids(num_angle_sets);

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
      // Check if the kernel is done
      if (streams_list[i].is_completed())
      {
        // Copy back outgoing (reflecting) boundary and non-local psi
        fluds_list[i]->CopyOutgoingPsiBackToHost(
          cbcd_sweep_chunk, angle_sets[i], in_flight_cell_ids[i]);
        // Update task dependencies
        auto& current_task_list = angle_sets[i]->GetCurrentTaskList();
        for (auto* task : in_flight_tasks[i])
        {
          for (uint64_t succ : task->successors)
          {
            --current_task_list[succ].num_dependencies;
            if (current_task_list[succ].num_dependencies == 0 and boundary_data_set[i])
              ready_queues[i].push_back(&current_task_list[succ]);
          }
          task->completed = true;
        }
        num_completed_tasks[i] += in_flight_tasks[i].size();
        // Send MPI data
        auto* comm = dynamic_cast<CBCD_AsynchronousCommunicator*>(angle_sets[i]->GetCommunicator());
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
      auto* comm = dynamic_cast<CBCD_AsynchronousCommunicator*>(angle_sets[i]->GetCommunicator());
      auto& current_task_list = angle_sets[i]->GetCurrentTaskList();
      auto received = comm->ReceiveData();
      if (not received.empty())
      {
        for (uint64_t t : received)
        {
          --current_task_list[t].num_dependencies;
          if (current_task_list[t].num_dependencies == 0 and boundary_data_set[i])
            ready_queues[i].push_back(&current_task_list[t]);
        }
        any_work_done = true;
      }
      comm->SendData();
    }

    // Set boundary data
    for (size_t i = 0; i < num_angle_sets; ++i)
    {
      if (executed[i] or boundary_data_set[i] or kernel_in_flight[i])
        continue;
      if (angle_sets[i]->IsDependencyResolved())
      {
        fluds_list[i]->CopyIncomingBoundaryPsiToDevice(cbcd_sweep_chunk, angle_sets[i]);
        boundary_data_set[i] = true;
        any_work_done = true;

        auto& current_task_list = angle_sets[i]->GetCurrentTaskList();
        for (auto& task : current_task_list)
        {
          if (task.num_dependencies == 0 and not task.completed)
            ready_queues[i].push_back(&task);
        }
      }
    }

    // Collect ready tasks and launch kernels (only if task dependencies changed)
    if (any_work_done)
    {
      for (size_t i = 0; i < num_angle_sets; ++i)
      {
        if (executed[i] or (not boundary_data_set[i]) or kernel_in_flight[i])
          continue;

        if (ready_queues[i].empty())
          continue;

        ready_tasks[i] = std::move(ready_queues[i]);
        ready_queues[i].clear();

        ready_cell_ids[i].clear();
        for (auto* task : ready_tasks[i])
          ready_cell_ids[i].push_back(task->reference_id);

        fluds_list[i]->CopyIncomingNonlocalPsiToDevice(angle_sets[i], ready_cell_ids[i]);
        cbcd_sweep_chunk.Sweep(ready_cell_ids[i], i);
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
      auto* comm = dynamic_cast<CBCD_AsynchronousCommunicator*>(angle_sets[i]->GetCommunicator());
      bool all_done = (num_completed_tasks[i] == current_task_list.size());
      if (all_done and comm->SendData())
      {
        angle_sets[i]->UpdateDependencyCounters();
        executed[i] = true;
        ++executed_anglesets;
        fluds_list[i]->CopySavedPsiFromDevice();
        auto* fluds = fluds_list[i];
        auto* as = angle_sets[i];
        // Cast away constness to add a callback
        streams_list[i].add_callback(
          [fluds, &cbcd_sweep_chunk, as]()
          { fluds->CopySavedPsiToDestinationPsi(cbcd_sweep_chunk, as); });
      }
    }
  }

  // Wait for the saved angular flux callbacks. CUDA and HIP order them before the copies below
  // through the legacy null stream, but SYCL queues are independent of each other.
  for (auto& stream : streams_list)
    stream.synchronize();

  // Copy phi and outflow data back to host
  cbcd_sweep_chunk.GetProblem().CopyPhiAndOutflowBackToHost();

  // Receive delayed data
  opensn::mpi_comm.barrier();
  bool received_delayed_data = false;
  while (not received_delayed_data)
  {
    received_delayed_data = true;

    for (const auto& angle_set : angle_sets)
    {
      if (angle_set->FlushSendBuffers() == AngleSetStatus::MESSAGES_PENDING)
        received_delayed_data = false;

      if (not angle_set->ReceiveDelayedData())
        received_delayed_data = false;
    }
  }

  // Reset all
  for (const auto& angle_set : angle_sets)
    angle_set->ResetSweepBuffers();
}

} // namespace opensn
