// SPDX-FileCopyrightText: 2026 The OpenSn Authors <https://open-sn.github.io/opensn/>
// SPDX-License-Identifier: MIT

#pragma once

#include <atomic>
#include <condition_variable>
#include <cstddef>
#include <functional>
#include <memory>
#include <mutex>
#include <semaphore>
#include <thread>
#include <vector>

namespace opensn
{

template <typename... Bases>
struct alignas(std::hardware_destructive_interference_size) CacheLineAligned : public Bases...
{
  template <typename... Args>
    requires(sizeof...(Args) == sizeof...(Bases)) &&
            (std::is_constructible_v<Bases, const Args&> && ...)
  constexpr explicit CacheLineAligned(const Args&... args) : Bases(args)...
  {
  }

  template <typename... Args>
    requires(sizeof...(Args) == sizeof...(Bases)) &&
            (std::is_constructible_v<Bases, Args &&> && ...)
  constexpr explicit CacheLineAligned(Args&&... args) : Bases(std::forward<Args>(args))...
  {
  }
};

/**
 * Single-Program Multiple-Data (SPMD) thread pool.
 * A fixed set of worker threads repeatedly executes the same callable in "epochs".
 * Each epoch publishes one task and wakes all workers. Each worker calls the task with its own rank
 * (thread index). The caller thread waits until all workers finish the current epoch.
 *
 * This design avoids repeated thread creation/join overhead when the number of workers is stable.
 *
 * \note Unless otherwise stated, public member functions are not thread-safe with respect to
 * concurrent calls from multiple external threads.
 */
class SPMD_ThreadPool
{
public:
  SPMD_ThreadPool() = default;

  /// Constructor with a fixed number of worker threads.
  explicit SPMD_ThreadPool(std::size_t n);

  SPMD_ThreadPool(const SPMD_ThreadPool&) = delete;
  SPMD_ThreadPool& operator=(const SPMD_ThreadPool&) = delete;

  ~SPMD_ThreadPool();

  std::size_t GetSize() const noexcept { return worker_threads_.size(); }

  /// Resize the numnber of worker threads.
  void Resize(std::size_t n);

  /**
   * Stop workers and join threads.
   * \note This method is safe to call multiple times.
   */
  void Stop();

  /**
   * Assign task to be run by all worker threads in the next epoch.
   * \warning Do not call this method while another thread is running. This method does not
   * guarantee thread safety.
   */
  template <class F>
  void AssignTask(F&& task)
  {
    std::scoped_lock<std::mutex> lock(mutex_);
    task_ = std::function<void(std::size_t)>(std::forward<F>(task));
  }
  /// Run the currently assigned task for the worker with the given index.
  void Run(std::size_t thread_idx);
  /// Wait until all workers have completed the current epoch.
  void WaitAll();

  /// Execute the currently assigned task in a new epoch.
  template <class F>
  void ExecuteBatch(F&& task)
  {
    AssignTask(std::forward<F>(task));
    const std::size_t n = worker_threads_.size();
    outstanding_.fetch_add((long long)n, std::memory_order_relaxed);
    for (std::size_t i = 0; i < n; ++i)
      signals_to_worker_[i]->release();
    WaitAll();
  }

private:
  /// Resize the worker vector and have each thread waiting on an infinite loop.
  void Start(std::size_t n);

  /// Run an infinite loop for the current calling until disrupted.
  void InfiniteLoop(std::size_t thread_idx);

  /// Persistent worker threads owned by the pool.
  std::vector<std::thread> worker_threads_;
  /// Per-worker semaphores for waking up workers to start an epoch.
  std::vector<std::unique_ptr<CacheLineAligned<std::binary_semaphore>>> signals_to_worker_;
  /// Mutex protecting shared state for task publication and epoch completion.
  std::mutex mutex_;
  /// Notifier that the current epoch has completed (i.e., all workers have finished).
  std::condition_variable cv_done_;

  /// Counter for the number of workers that have not yet completed the current epoch.
  std::atomic<long long> outstanding_{0};
  /// Flag indicating whether the thread pool is stopping.
  std::atomic<bool> stopped_{false};

  /**
   * Published per-epoch task.
   * \note This callable is valid only while an epoch is in progress.
   */
  std::function<void(std::size_t)> task_;
};

} // namespace opensn
