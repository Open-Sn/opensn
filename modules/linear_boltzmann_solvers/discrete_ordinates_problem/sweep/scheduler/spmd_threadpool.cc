// SPDX-FileCopyrightText: 2026 The OpenSn Authors <https://open-sn.github.io/opensn/>
// SPDX-License-Identifier: MIT

#include "modules/linear_boltzmann_solvers/discrete_ordinates_problem/sweep/scheduler/spmd_threadpool.h"
#include <cassert>

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
  if (worker_threads_.size() == n)
    return;
  Stop();
  if (n > 0)
    Start(n);
}

void
SPMD_ThreadPool::Start(std::size_t n)
{
  stopped_.store(false, std::memory_order_release);
  signals_to_worker_.reserve(n);
  worker_threads_.reserve(n);
  for (std::size_t i = 0; i < n; ++i)
    signals_to_worker_.emplace_back(std::make_unique<CacheLineAligned<std::binary_semaphore>>(0));
  for (std::size_t i = 0; i < n; ++i)
    worker_threads_.emplace_back([this, i] { InfiniteLoop(i); });
}

void
SPMD_ThreadPool::Stop()
{
  bool expected = false;
  if (!stopped_.compare_exchange_strong(expected, true, std::memory_order_acq_rel))
    return;

  for (auto& s : signals_to_worker_)
    s->release();

  for (auto& t : worker_threads_)
    if (t.joinable())
      t.join();

  worker_threads_.clear();
  signals_to_worker_.clear();
  task_ = nullptr;
  outstanding_.store(0, std::memory_order_release);
}

void
SPMD_ThreadPool::Run(std::size_t thread_idx)
{
  assert(thread_idx < worker_threads_.size());
  assert(task_ != nullptr);
  outstanding_.fetch_add(1, std::memory_order_relaxed);
  signals_to_worker_[thread_idx]->release();
}

void
SPMD_ThreadPool::WaitAll()
{
  std::unique_lock<std::mutex> lock(mutex_);
  cv_done_.wait(lock, [this] { return outstanding_.load(std::memory_order_acquire) == 0; });
}

void
SPMD_ThreadPool::InfiniteLoop(std::size_t thread_idx)
{
  while (true)
  {
    signals_to_worker_[thread_idx]->acquire();
    if (stopped_.load(std::memory_order_acquire))
      return;

    if (task_ != nullptr)
      task_(thread_idx);

    if (outstanding_.fetch_sub(1, std::memory_order_acq_rel) == 1)
      cv_done_.notify_all();
  }
}

} // namespace opensn
