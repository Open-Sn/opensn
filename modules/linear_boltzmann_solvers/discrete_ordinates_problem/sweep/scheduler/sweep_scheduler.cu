// SPDX-FileCopyrightText: 2026 The OpenSn Authors <https://open-sn.github.io/opensn/>
// SPDX-License-Identifier: MIT

#include "modules/linear_boltzmann_solvers/discrete_ordinates_problem/sweep/scheduler/sweep_scheduler.h"
#include "modules/linear_boltzmann_solvers/discrete_ordinates_problem/sweep/angle_set/aahd_angle_set.h"
#include "modules/linear_boltzmann_solvers/discrete_ordinates_problem/sweep_chunks/aahd_sweep_chunk.h"
#include "modules/linear_boltzmann_solvers/discrete_ordinates_problem/sweep_chunks/cbcd_sweep_chunk.h"
#include "modules/linear_boltzmann_solvers/discrete_ordinates_problem/discrete_ordinates_problem.h"
#include "caribou/main.hpp"
#include "caliper/cali.h"
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
  CALI_CXX_MARK_SCOPE("SweepScheduler::ScheduleAlgoAsyncFIFO");

  auto& cbcd_sweep_chunk = dynamic_cast<CBCDSweepChunk&>(sweep_chunk);
  cbcd_sweep_chunk.GetProblem().CopyPhiAndSrcToDevice();
  const auto groupset_id = angle_agg_.GetGroupsetID();
  if (cbcd_sweep_chunk.NeedsBoundaryUpload())
    cbcd_sweep_chunk.GetProblem().TransferDeviceBoundaryData(groupset_id, true);
  cbcd_sweep_chunk.RefreshKernelArguments();

  const auto& angle_sets = cbcd_sweep_chunk.GetAngleSets();
  const auto num_angle_sets = angle_sets.size();
  for (auto* angle_set : angle_sets)
    angle_set->ResetSweepDependencies();

  // Worker w sweeps the angle sets that CBCD_AngleSetWorker assigns to it. The communicator sizes
  // its worker queues for the same assignment. Each worker visit advances each of its
  // unfinished angle sets by one step and then launches the cells they made ready in one fused
  // launch.
  const auto num_workers = pool_.GetSize();
  cbcd_sweep_chunk.StartCommunicator(num_workers);
  pool_.ExecuteBatch(
    [num_workers, num_angle_sets, &angle_sets, &cbcd_sweep_chunk](std::size_t worker_id)
    {
      std::vector<std::size_t> active_angle_set_ids;
      active_angle_set_ids.reserve((num_angle_sets + num_workers - 1) / num_workers);
      for (std::size_t i = 0; i < num_angle_sets; ++i)
        if (CBCD_AngleSetWorker(i, num_workers) == worker_id)
          active_angle_set_ids.push_back(i);

      while (not active_angle_set_ids.empty())
      {
        cbcd_sweep_chunk.PollLane(worker_id);
        bool any_work_done = false;
        for (std::size_t i = 0; i < active_angle_set_ids.size();)
        {
          auto* angle_set = angle_sets[active_angle_set_ids[i]];
          if (not angle_set->IsSweepInitialized())
          {
            if (not angle_set->TryInitialize())
            {
              ++i;
              continue;
            }
            any_work_done = true;
          }

          any_work_done |= angle_set->TryAdvanceOneStep(cbcd_sweep_chunk, worker_id);

          if (angle_set->IsSweepComplete())
          {
            active_angle_set_ids[i] = active_angle_set_ids.back();
            active_angle_set_ids.pop_back();
            continue;
          }

          ++i;
        }
        cbcd_sweep_chunk.SubmitFusedLaunch(worker_id);
        if (not any_work_done)
          std::this_thread::yield();
      }
    });

  cbcd_sweep_chunk.StopCommunicator();
  opensn::mpi_comm.barrier();

  cbcd_sweep_chunk.GetProblem().TransferDeviceBoundaryData(groupset_id, false);
  cbcd_sweep_chunk.GetProblem().CopyPhiAndOutflowBackToHost();

  for (auto* angle_set : angle_sets)
    angle_set->ResetSweepBuffers();
}

} // namespace opensn
