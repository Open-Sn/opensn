// SPDX-FileCopyrightText: 2026 The OpenSn Authors <https://open-sn.github.io/opensn/>
// SPDX-License-Identifier: MIT

#pragma once

#include "modules/linear_boltzmann_solvers/discrete_ordinates_problem/sweep/scheduler/spmd_threadpool.h"
#include <atomic>
#include <cassert>
#include <cstddef>
#include <span>
#include <vector>

namespace opensn
{

static_assert(HardwareInterferenceSize != 0 and
                (HardwareInterferenceSize & (HardwareInterferenceSize - 1)) == 0,
              "Device CBC queue alignment must be a nonzero power of two.");

/// Append-only host queue with one producer and one consumer per sweep.
/// Size it for all elements in a sweep. Entries are neither overwritten nor recycled before Reset.
template <typename T>
// Deliberate padding keeps producer and consumer indices on separate cache lines.
class LockFreeSPSCQueue // NOLINT(clang-analyzer-optin.performance.Padding)
{
public:
  /// Allocate storage for capacity elements between resets.
  explicit LockFreeSPSCQueue(std::size_t capacity) : elements_(capacity) {}

  /// Remove all elements. Call only while neither end is in use.
  void Reset()
  {
    published_.store(0, std::memory_order_relaxed);
    consumed_ = 0;
  }

  /// Append from the sole producer with capacity remaining. Release-publish preceding writes.
  void Push(const T& element)
  {
    const auto size = published_.load(std::memory_order_relaxed);
    assert(size < elements_.size());
    elements_[size] = element;
    published_.store(size + 1, std::memory_order_release);
  }

  /// Acquire entries published since the sole consumer's previous call, with preceding producer
  /// writes. The returned span borrows storage and remains valid until Reset or destruction.
  std::span<const T> Take()
  {
    const auto size = published_.load(std::memory_order_acquire);
    const std::span<const T> ready(elements_.data() + consumed_, size - consumed_);
    consumed_ = size;
    return ready;
  }

  /// Check for unread entries at this snapshot. Only the consumer may call this.
  bool IsEmpty() const { return published_.load(std::memory_order_acquire) == consumed_; }

private:
  std::vector<T> elements_;
  alignas(HardwareInterferenceSize) std::atomic<std::size_t> published_{0};
  alignas(HardwareInterferenceSize) std::size_t consumed_ = 0;
};

} // namespace opensn
