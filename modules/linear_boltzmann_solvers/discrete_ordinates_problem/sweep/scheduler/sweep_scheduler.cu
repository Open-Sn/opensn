// SPDX-FileCopyrightText: 2026 The OpenSn Authors <https://open-sn.github.io/opensn/>
// SPDX-License-Identifier: MIT

#include "modules/linear_boltzmann_solvers/discrete_ordinates_problem/sweep/scheduler/sweep_scheduler.h"
#include "modules/linear_boltzmann_solvers/discrete_ordinates_problem/sweep/angle_set/aahd_angle_set.h"
#include "modules/linear_boltzmann_solvers/discrete_ordinates_problem/sweep_chunks/aahd_sweep_chunk.h"
#include "modules/linear_boltzmann_solvers/discrete_ordinates_problem/discrete_ordinates_problem.h"
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

  // reset dependency counters and pre-post receives for all anglesets
  std::size_t num_anglesets = angle_agg_.GetNumAngleSets();
  execution_order_.resize(num_anglesets);
  for (std::size_t i = 0; i < num_anglesets; ++i)
  {
    auto aahd_angle_set = static_cast<AAHD_AngleSet*>(angle_agg_[i].get());
    aahd_angle_set->ResetDependencyCounter();
    aahd_angle_set->PrepostReceives();
    execution_order_[i] = i;
  }

  // assign sweep task to thread pool (but not execution yet)
  pool_.AssignTask(
    [this, &sweep_chunk](std::size_t i)
    {
      auto* aahd = static_cast<AAHD_AngleSet*>(angle_agg_[i].get());
      aahd->AngleSetAdvance(sweep_chunk, AngleSetStatus::EXECUTE);
    });

  // poll for readiness and launch threads
  while (!execution_order_.empty())
  {
    for (auto it = execution_order_.begin(); it != execution_order_.end();)
    {
      auto* angle_set = static_cast<AAHD_AngleSet*>(angle_agg_[*it].get());
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
    auto aahd_angle_set = static_cast<AAHD_AngleSet*>(angle_set.get());
    aahd_angle_set->WaitForDownstreamAndDelayed();
  }

  // copy phi and outflow data back to host
  aah_sweep_chunk.GetProblem().CopyPhiAndOutflowBackToHost();

  // reset all anglesets
  for (auto& angle_set : angle_agg_)
    angle_set->ResetSweepBuffers();
}

} // namespace opensn
