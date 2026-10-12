// SPDX-FileCopyrightText: 2026 The OpenSn Authors <https://open-sn.github.io/opensn/>
// SPDX-License-Identifier: MIT

#include "modules/linear_boltzmann_solvers/discrete_ordinates_problem/sweep/scheduler/spmd_threadpool.h"
#include "framework/runtime.h"
#include <cassert>

namespace opensn
{
namespace
{

std::size_t
ResolveWorkerCount(std::size_t requested_workers)
{
  if (requested_workers == 0)
    return 0;

  // Only cap when OPENSN_NUM_THREADS is explicitly set; otherwise use the
  // requested count so callers that size the pool to their work item count
  // are not silently throttled to one thread.
  if (std::getenv("OPENSN_NUM_THREADS") == nullptr) // NOLINT(concurrency-mt-unsafe)
    return requested_workers;

  return std::min<std::size_t>(requested_workers, std::max(1U, opensn_num_threads));
}

} // namespace

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
  n = ResolveWorkerCount(n);
  if (worker_threads_.size() == n)
    return;
  Stop();
  if (n > 0)
    Start(n);
}

void
SPMD_ThreadPool::Start(std::size_t n)
{
  if (workers_initialized_ || n == 0)
    return;

  worker_threads_.reserve(n);
  epoch_states_.assign(n, EpochState{0, 0});
  outstanding_ = 0;

  for (std::size_t i = 0; i < n; ++i)
    worker_threads_.emplace_back(&SPMD_ThreadPool::InfiniteLoop, this, i);

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
  }
  cv_start_.notify_all();
  cv_done_.notify_all();

  for (auto& t : worker_threads_)
    if (t.joinable())
      t.join();

  worker_threads_.clear();
  epoch_states_.clear();
  task_ = nullptr;

  workers_initialized_ = false;
  stop_workers_ = false;
  outstanding_ = 0;
}

void
SPMD_ThreadPool::Run(std::size_t thread_idx)
{
  assert(workers_initialized_);
  assert(thread_idx < worker_threads_.size());
  {
    std::scoped_lock<std::mutex> lock(mutex_);
    ++epoch_states_[thread_idx].request;
    ++outstanding_;
  }
  cv_start_.notify_all();
}

void
SPMD_ThreadPool::WaitAll()
{
  std::unique_lock<std::mutex> lock(mutex_);
  cv_done_.wait(lock, [this] { return outstanding_ == 0; });
}

void
SPMD_ThreadPool::InfiniteLoop(std::size_t thread_idx)
{
  while (true)
  {
    {
      std::unique_lock<std::mutex> lock(mutex_);
      cv_start_.wait(lock,
                     [&] { return stop_workers_ || epoch_states_[thread_idx].IsIncomplete(); });
      if (stop_workers_)
        return;
      ++epoch_states_[thread_idx].done;
    }
    task_(thread_idx);
    {
      std::scoped_lock<std::mutex> lock(mutex_);
      if (--outstanding_ == 0)
        cv_done_.notify_all();
    }
  }
}

} // namespace opensn
