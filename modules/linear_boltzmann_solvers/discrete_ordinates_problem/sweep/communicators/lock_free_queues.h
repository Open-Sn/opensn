// SPDX-FileCopyrightText: 2026 The OpenSn Authors <https://open-sn.github.io/opensn/>
// SPDX-License-Identifier: MIT

#pragma once

#include <algorithm>
#include <atomic>
#include <cstddef>
#include <limits>
#include <new>
#include <thread>
#include <vector>

namespace opensn
{

/** Bounded lock-free SPSC ring whose elements retain their allocated capacity. */
template <typename T>
class LockFreeSPSCSlotQueue
{
public:
  /** Allocates the ring storage and resets both endpoints. */
  void Preallocate(const std::size_t capacity)
  {
    buffer_ = std::vector<T>(capacity);
    producer_head_ = 0;
    consumer_tail_ = 0;
    producer_tail_cache_ = 0;
    consumer_head_cache_ = 0;
    published_head_.store(0, std::memory_order_relaxed);
    consumed_tail_.store(0, std::memory_order_relaxed);
  }

  /** Applies `callback` once to each reusable element. */
  template <typename Callback>
  void InitializeSlots(Callback&& callback)
  {
    for (auto& element : buffer_)
      callback(element);
  }

  /** Returns the next producer-owned element, waiting only while the ring is full. */
  T& ReserveSlot()
  {
    const auto capacity = buffer_.size();
    while ((producer_head_ - producer_tail_cache_) >= capacity)
    {
      producer_tail_cache_ = consumed_tail_.load(std::memory_order_acquire);
      if ((producer_head_ - producer_tail_cache_) < capacity)
        break;
      std::this_thread::yield();
    }
    return buffer_[producer_head_ % capacity];
  }

  /** Makes the most recently reserved element visible to the consumer. */
  void PublishSlot()
  {
    ++producer_head_;
    published_head_.store(producer_head_, std::memory_order_release);
  }

  /** Returns pointers to the current contiguous logical prefix of ready elements. */
  void PeekReadySlots(std::vector<T*>& out,
                      const std::size_t max_slots = std::numeric_limits<std::size_t>::max())
  {
    out.clear();
    if (buffer_.empty())
      return;

    consumer_head_cache_ = published_head_.load(std::memory_order_acquire);
    const auto ready_count = std::min(consumer_head_cache_ - consumer_tail_, max_slots);
    out.reserve(ready_count);
    const auto capacity = buffer_.size();
    for (std::size_t i = 0; i < ready_count; ++i)
      out.push_back(&buffer_[(consumer_tail_ + i) % capacity]);
  }

  /** Releases `count` elements previously returned by `PeekReadySlots`. */
  void ReleaseReadySlots(const std::size_t count)
  {
    if (count == 0)
      return;

    consumer_tail_ += count;
    consumed_tail_.store(consumer_tail_, std::memory_order_release);
  }

  /** Processes and releases every element currently visible to the consumer. */
  template <typename Callback>
  std::size_t ProcessReady(Callback&& callback)
  {
    if (buffer_.empty())
      return 0;

    consumer_head_cache_ = published_head_.load(std::memory_order_acquire);
    const auto ready_count = consumer_head_cache_ - consumer_tail_;
    const auto capacity = buffer_.size();
    for (std::size_t i = 0; i < ready_count; ++i)
      callback(buffer_[(consumer_tail_ + i) % capacity]);
    ReleaseReadySlots(ready_count);
    return ready_count;
  }

  /** Returns whether the consumer has no published elements. */
  bool Empty() const
  {
    if (buffer_.empty())
      return true;
    return published_head_.load(std::memory_order_acquire) == consumer_tail_;
  }

private:
  std::vector<T> buffer_;
  alignas(std::hardware_destructive_interference_size) std::size_t producer_head_ = 0;
  alignas(std::hardware_destructive_interference_size) std::size_t consumer_tail_ = 0;
  std::size_t producer_tail_cache_ = 0;
  std::size_t consumer_head_cache_ = 0;
  alignas(std::hardware_destructive_interference_size) std::atomic<std::size_t> published_head_{0};
  alignas(std::hardware_destructive_interference_size) std::atomic<std::size_t> consumed_tail_{0};
};

} // namespace opensn
