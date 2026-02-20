// SPDX-FileCopyrightText: 2026 The OpenSn Authors <https://open-sn.github.io/opensn/>
// SPDX-License-Identifier: MIT

#pragma once

#include <condition_variable>
#include <functional>
#include <cstddef>
#include <exception>
#include <mutex>
#include <thread>
#include <vector>

namespace opensn
{

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
  bool IsInitialized() const noexcept { return workers_initialized_; }

  /// Resize the numnber of worker threads.
  void Resize(std::size_t n);

  /**
   * Stop workers and join threads.
   * \note This method is safe to call multiple times.
   */
  void Stop();

  /// Run one SPMD epoch, then wait for all tasks to finish.
  template <class F>
  void run(F&& task)
  {
    if (!workers_initialized_)
      Start(0);
    if (worker_threads_.empty())
      throw std::runtime_error("SPMD_ThreadPool has size 0");

    // ensure that task stay alive after all workers have finished
    {
      std::scoped_lock<std::mutex> lock(mutex_);
      task_fn_ = std::function<void(std::size_t)>(std::forward<F>(task));
      workers_remaining_ = worker_threads_.size();
      ++work_epoch_;
    }
    cv_start_.notify_all();

    // wait for all workers to complete this epoch
    {
      std::unique_lock<std::mutex> lock(mutex_);
      cv_done_.wait(lock, [this] { return workers_remaining_ == 0; });
      task_fn_ = nullptr;
    }
  }

private:
  /// Resize the worker vector and have each thread waiting on an infinite loop.
  void Start(std::size_t n);

  /// Run an infinite loop for the current calling until disrupted.
  void InfiniteLoop(std::size_t thread_idx);

  /// Persistent worker threads owned by the pool.
  std::vector<std::thread> worker_threads_;
  /// Mutex protecting shared state for task publication and epoch completion.
  std::mutex mutex_;
  /// Notifier that a new epoch has started or that shutdown is requested.
  std::condition_variable cv_start_;
  /// Notifier that the current epoch has completed (i.e., all workers have finished).
  std::condition_variable cv_done_;

  /// Flag indicating whether workers have been started at least once.
  bool workers_initialized_ = false;
  /// Flag indicating whether workers should stop.
  bool stop_workers_ = false;
  /// Number of workers remaining to finish the current epoch.
  std::size_t workers_remaining_ = 0;
  /// Monotonic epoch counter used to distinguish work submissions and to avoid missed wakeups.
  std::size_t work_epoch_ = 0;

  /**
   * Published per-epoch task.
   * \note This callable is valid only while an epoch is in progress.
   */
  std::function<void(std::size_t)> task_fn_;
};

} // namespace opensn
