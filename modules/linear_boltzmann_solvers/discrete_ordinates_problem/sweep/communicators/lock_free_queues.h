// SPDX-FileCopyrightText: 2026 The OpenSn Authors <https://open-sn.github.io/opensn/>
// SPDX-License-Identifier: MIT

#pragma once

#include "modules/linear_boltzmann_solvers/discrete_ordinates_problem/sweep/scheduler/spmd_threadpool.h"
#include <algorithm>
#include <array>
#include <atomic>
#include <cassert>
#include <cstddef>
#include <span>
#include <stdexcept>
#include <thread>
#include <vector>

namespace opensn
{

/**
 * Bounded lock-free SPSC ring whose elements retain their allocated capacity.
 */
template <typename T>
class LockFreeSPSCSlotQueue
{
public:
  /// Allocate the ring storage and reset both endpoints.
  void Preallocate(const std::size_t capacity)
  {
    buffer_ = std::vector<T>(capacity);
    producer_ = {};
    producer_.data = buffer_.data();
    producer_.capacity = capacity;
    consumer_ = {};
    consumer_.data = buffer_.data();
    consumer_.capacity = capacity;
    published_head_.store(0, std::memory_order_relaxed);
    consumed_tail_.store(0, std::memory_order_relaxed);
  }

  /**
   * Return the next producer-owned element, which waits only while the ring is full.
   * Throws std::logic_error if the queue has no storage since it could never accept an element.
   */
  T& ReserveSlot()
  {
    const auto capacity = producer_.capacity;
    while ((producer_.head - producer_.tail_cache) >= capacity)
    {
      if (capacity == 0)
        throw std::logic_error("LockFreeSPSCSlotQueue: an element was reserved in a queue without "
                               "storage.");
      producer_.tail_cache = consumed_tail_.load(std::memory_order_acquire);
      if ((producer_.head - producer_.tail_cache) < capacity)
        break;
      std::this_thread::yield();
    }
    return producer_.data[producer_.index];
  }

  /// Make the most recently reserved element visible to the consumer.
  void PublishSlot()
  {
    if (++producer_.index == producer_.capacity)
      producer_.index = 0;
    ++producer_.head;
    published_head_.store(producer_.head, std::memory_order_release);
  }

  /// Return the published prefix as at most two spans, which is valid until released.
  std::array<std::span<const T>, 2> PeekReadySlots()
  {
    const auto capacity = consumer_.capacity;
    if (capacity == 0)
      return {};

    const auto ready_count = published_head_.load(std::memory_order_acquire) - consumer_.tail;
    if (ready_count == 0)
      return {};

    const auto first_index = consumer_.index;
    const auto first_count = std::min(ready_count, capacity - first_index);
    const std::span<const T> buffer(consumer_.data, capacity);
    return {buffer.subspan(first_index, first_count), buffer.first(ready_count - first_count)};
  }

  /// Release `count` elements previously returned by `PeekReadySlots`.
  void ReleaseReadySlots(const std::size_t count)
  {
    if (count == 0)
      return;

    assert(count <= consumer_.capacity);
    const auto until_wrap = consumer_.capacity - consumer_.index;
    consumer_.index = count >= until_wrap ? count - until_wrap : consumer_.index + count;
    consumer_.tail += count;
    consumed_tail_.store(consumer_.tail, std::memory_order_release);
  }

  /// Process and release every element currently visible to the consumer.
  template <typename Callback>
  std::size_t ProcessReady(Callback callback)
  {
    const auto ready = PeekReadySlots();
    const auto ready_count = ready[0].size() + ready[1].size();
    for (const auto span : ready)
      for (const auto& element : span)
        callback(element);
    ReleaseReadySlots(ready_count);
    return ready_count;
  }

  /// Return whether the consumer has no published elements.
  bool Empty() const
  {
    if (consumer_.capacity == 0)
      return true;
    return published_head_.load(std::memory_order_acquire) == consumer_.tail;
  }

private:
  struct alignas(HardwareInterferenceSize) ProducerState
  {
    T* data = nullptr;
    std::size_t capacity = 0;
    std::size_t head = 0;
    std::size_t tail_cache = 0;
    std::size_t index = 0;
  };
  struct alignas(HardwareInterferenceSize) ConsumerState
  {
    T* data = nullptr;
    std::size_t capacity = 0;
    std::size_t tail = 0;
    std::size_t index = 0;
  };

  ProducerState producer_;
  ConsumerState consumer_;
  alignas(HardwareInterferenceSize) std::atomic<std::size_t> published_head_{0};
  alignas(HardwareInterferenceSize) std::atomic<std::size_t> consumed_tail_{0};
  /// Element storage, which is used only through the endpoints' copies of its address.
  std::vector<T> buffer_;
};

} // namespace opensn
