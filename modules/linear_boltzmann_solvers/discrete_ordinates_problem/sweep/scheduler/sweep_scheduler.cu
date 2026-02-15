// SPDX-FileCopyrightText: 2026 The OpenSn Authors <https://open-sn.github.io/opensn/>
// SPDX-License-Identifier: MIT

#include "modules/linear_boltzmann_solvers/discrete_ordinates_problem/sweep/scheduler/sweep_scheduler.h"
#include "modules/linear_boltzmann_solvers/discrete_ordinates_problem/sweep/angle_set/aahd_angle_set.h"
#include "modules/linear_boltzmann_solvers/discrete_ordinates_problem/sweep_chunks/aahd_sweep_chunk.h"
#include "modules/linear_boltzmann_solvers/discrete_ordinates_problem/discrete_ordinates_problem.h"
#include "caliper/cali.h"

namespace opensn
{

void
SweepScheduler::ScheduleAlgoAAO(SweepChunk& sweep_chunk)
{
  CALI_CXX_MARK_SCOPE("SweepScheduler::ScheduleAlgoAAO");

  // copy phi and src moments to device
  auto aah_sweep_chunk = static_cast<AAHDSweepChunk&>(sweep_chunk);
  aah_sweep_chunk.GetProblem().CopyPhiAndSrcToDevice();

  std::size_t num_anglesets = angle_agg_.GetNumAngleSets();

  // Initialize (or reinitialize) persistent workers if needed.
  if (aao_workers_initialized_ && aao_worker_threads_.size() != num_anglesets)
  {
    {
      std::scoped_lock<std::mutex> lock(aao_mutex_);
      aao_stop_workers_ = true;
      ++aao_work_epoch_;
    }
    aao_cv_start_.notify_all();
    for (auto& thread : aao_worker_threads_)
      if (thread.joinable())
        thread.join();
    aao_worker_threads_.clear();
    aao_workers_initialized_ = false;
    aao_stop_workers_ = false;
    aao_workers_remaining_ = 0;
    aao_active_sweep_chunk_ = nullptr;
    aao_work_epoch_ = 0;
  }

  if (not aao_workers_initialized_)
  {
    aao_worker_threads_.reserve(num_anglesets);
    for (std::size_t i = 0; i < num_anglesets; ++i)
    {
      aao_worker_threads_.emplace_back(
        [this, i]()
        {
          auto aahd_angle_set = static_cast<AAHD_AngleSet*>(angle_agg_[i].get());
          std::size_t seen_epoch = 0;
          while (true)
          {
            SweepChunk* active_chunk = nullptr;
            {
              std::unique_lock<std::mutex> lock(aao_mutex_);
              aao_cv_start_.wait(lock,
                                 [this, seen_epoch]()
                                 { return aao_stop_workers_ or aao_work_epoch_ > seen_epoch; });
              if (aao_stop_workers_)
                return;
              seen_epoch = aao_work_epoch_;
              active_chunk = aao_active_sweep_chunk_;
            }

            aahd_angle_set->AngleSetAdvance(*active_chunk, AngleSetStatus::EXECUTE);

            {
              std::scoped_lock<std::mutex> lock(aao_mutex_);
              if (aao_workers_remaining_ > 0)
                --aao_workers_remaining_;
              if (aao_workers_remaining_ == 0)
                aao_cv_done_.notify_one();
            }
          }
        });
    }
    aao_workers_initialized_ = true;
  }

  // set the latches for all anglesets
  for (auto& angle_set : angle_agg_)
  {
    auto aahd_angle_set = static_cast<AAHD_AngleSet*>(angle_set.get());
    aahd_angle_set->SetStartingLatch();
  }

  {
    std::scoped_lock<std::mutex> lock(aao_mutex_);
    aao_active_sweep_chunk_ = &sweep_chunk;
    aao_workers_remaining_ = num_anglesets;
    ++aao_work_epoch_;
  }
  aao_cv_start_.notify_all();

  // wait for all workers to complete this sweep
  {
    std::unique_lock<std::mutex> lock(aao_mutex_);
    aao_cv_done_.wait(lock, [this]() { return aao_workers_remaining_ == 0; });
  }

  // copy phi and outflow data back to host
  aah_sweep_chunk.GetProblem().CopyPhiAndOutflowBackToHost();

  // reset all anglesets
  for (auto& angle_set : angle_agg_)
    angle_set->ResetSweepBuffers();
}

} // namespace opensn
