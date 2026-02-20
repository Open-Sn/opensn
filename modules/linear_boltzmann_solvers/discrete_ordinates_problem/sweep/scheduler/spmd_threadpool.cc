// SPDX-FileCopyrightText: 2026 The OpenSn Authors <https://open-sn.github.io/opensn/>
// SPDX-License-Identifier: MIT

#include "modules/linear_boltzmann_solvers/discrete_ordinates_problem/sweep/scheduler/spmd_threadpool.h"

namespace opensn
{

SPMD_ThreadPool::SPMD_ThreadPool(std::size_t n)
{
  Resize(n);
}

SPMD_ThreadPool::~SPMD_ThreadPool()
{
  Stop();
}

void
SPMD_ThreadPool::Resize(std::size_t n)
{
  if (workers_initialized_ && worker_threads_.size() != n)
    Stop();
  if (!workers_initialized_ && n > 0)
    Start(n);
}

void
SPMD_ThreadPool::Start(std::size_t n)
{
  if (workers_initialized_)
    return;
  if (n == 0)
    return;

  worker_threads_.reserve(n);
  for (std::size_t thread_idx = 0; thread_idx < n; ++thread_idx)
  {
    worker_threads_.emplace_back(&SPMD_ThreadPool::InfiniteLoop, this, thread_idx);
  }
  workers_initialized_ = true;
}

void
SPMD_ThreadPool::Stop()
{
  if (!workers_initialized_)
    return;

  {
    std::scoped_lock<std::mutex> lock(mutex_);
    stop_workers_ = true;
    ++work_epoch_; // ensure wait predicate flips
  }
  cv_start_.notify_all();

  for (auto& t : worker_threads_)
    if (t.joinable())
      t.join();

  worker_threads_.clear();

  workers_initialized_ = false;
  stop_workers_ = false;
  workers_remaining_ = 0;
  work_epoch_ = 0;
  task_fn_ = nullptr;
}

void
SPMD_ThreadPool::InfiniteLoop(std::size_t thread_idx)
{
  std::size_t seen_epoch = 0;
  while (true)
  {
    std::function<void(std::size_t)>* local_fn = nullptr;
    {
      std::unique_lock<std::mutex> lock(mutex_);
      cv_start_.wait(lock,
                     [this, &seen_epoch] { return stop_workers_ || work_epoch_ > seen_epoch; });

      if (stop_workers_)
        return;

      seen_epoch = work_epoch_;
      local_fn = &task_fn_;
    }
    (*local_fn)(thread_idx);
    {
      std::scoped_lock<std::mutex> lock(mutex_);
      if (workers_remaining_ > 0)
        --workers_remaining_;
      if (workers_remaining_ == 0)
        cv_done_.notify_one();
    }
  }
}

} // namespace opensn
